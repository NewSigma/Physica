/*
 * Copyright 2023-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/CuboidLinear.h"
#include "PeriodicCell.h"
#include "PeriodIndex3D.h"

namespace Physica {
    /**
     * References:
     * [1] Phys. Rev. 178 (1969) 1419–1423; http://dx.doi.org/10.1103/PhysRev.178.1419.
     */
    template<Scalar T>
    class BrillouInterpolate {
        constexpr static unsigned int Dim = 3;
        constexpr static bool isComplex = T::isComplex();
        using This = BrillouInterpolate<T>;
        using Tr = T::RealType;
        using Tc = T::ComplexType;
        using LatticeMatrix = PeriodicCell<Tr, Dim>::LatticeMatrix;
        using CoeffMatrixM = DenseMatrix<Tc, MatrixMajor::Row>;
        using Vec3D = Vector3D<Tr>;

        DenseTensor<Tc, 3> baseCoeff;
        LatticeMatrix lattice;
        DenseLU<Tc, true> lu;
        VectorND<Tc> solverBuffer;
        Tr smoothFactor1;
        Tr smoothFactor2;
        Index3D dataDim;
    public:
        BrillouInterpolate(Index3D baseDim, LatticeMatrix unitcell, Tr smoothFactor1_, Tr smoothFactor2_);
        BrillouInterpolate(const This& obj) = default;
        BrillouInterpolate(This&& obj) noexcept = default;
        ~BrillouInterpolate() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T operator()(Vec3D kPoint) const;
        /* Operations */
        Tr calcRoughness() const;
        void interpolate(const DenseTensor<T, 3>& data);
        T interpolateFEM(Vec3D kPoint, const DenseTensor<T, 3>& data) const;
        template<Matrix M>
        M interpolateHermiteMatrix(Vec3D qPoint, const ArrayND<M, 3>& matrixGrid);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Index3D getBaseDim() const noexcept { return baseCoeff.getShape(); }
    private:
        void initBaseCoeff();
        CoeffMatrixM makeMatrixM() const;
    };

    template<Scalar T>
    BrillouInterpolate<T>::BrillouInterpolate(Index3D baseDim, LatticeMatrix unitcell, Tr smoothFactor1_, Tr smoothFactor2_)
            : baseCoeff(baseDim)
            , lattice(unitcell)
            , smoothFactor1(smoothFactor1_)
            , smoothFactor2(smoothFactor2_)
            , dataDim({0, 0, 0}) {
        Tr minSquaredNorm = std::numeric_limits<Tr>::max();
        for (int i = 0; i < 3; ++i)
            minSquaredNorm = std::min(minSquaredNorm, lattice.row(i).squaredNorm());
        smoothFactor1 /= minSquaredNorm;
        smoothFactor2 /= square(minSquaredNorm);
    }

    template<Scalar T>
    T BrillouInterpolate<T>::operator()(Vec3D kPoint) const {
        T result(0);
        baseCoeff.forND([this, kPoint, &result](const auto& coeff, Index3D index) {
            Tr phase(0);
            for (size_t dim = 0; dim < Dim; ++dim)
                phase += Tr(index[dim]) * kPoint[dim];
            phase *= Tr(2 * M_PI);
            const auto factor = cos(phase);

            if constexpr (!isComplex) {
                PeriodIndex3D pIndex(index, getBaseDim());
                if (pIndex.isInReducedK())
                    result += (coeff * factor).real();
                else
                    result += (baseCoeff[pIndex.toReducedK()].conjugate() * factor).real();
            }
            else {
                (void)this; // Silence unused 'this' warning;
                result += coeff * factor;
            }
        });
        return result;
    }

    template<Scalar T>
    auto BrillouInterpolate<T>::calcRoughness() const -> Tr {
        Tr result = 0;
        const auto kernel = [this, &result](Vec3D r, Index3D index) {
            const Tr r2 = r.squaredNorm();
            const bool isGammaPoint = r2 < std::numeric_limits<Tr>::epsilon();
            if (isGammaPoint)
                return;
            const Tr r4 = square(r2);
            result += baseCoeff[index].squaredNorm() * (Tr(1) + smoothFactor1 * r2 + smoothFactor2 * r4);
        };
        PeriodicCell<Tr, 3>::template forLatticeCloud<true>(kernel, lattice, getBaseDim());
        return result;
    }

    template<Scalar T>
    void BrillouInterpolate<T>::interpolate(const DenseTensor<T, 3>& data) {
        assert(dataDim[0] <= getBaseDim()[0] && "[Error]: Not enough base function, interpolation is overdetermined");
        assert(dataDim[1] <= getBaseDim()[1] && "[Error]: Not enough base function, interpolation is overdetermined");
        assert(dataDim[2] <= getBaseDim()[2] && "[Error]: Not enough base function, interpolation is overdetermined");
        const bool useNormalFFT = dataDim[0] == getBaseDim()[0] && dataDim[1] == getBaseDim()[1] && dataDim[2] == getBaseDim()[2];
        if (useNormalFFT) {
            dataDim = data.getShape();
            FFT<Tc, 3> fft(dataDim, PlanFlag::Estimate);
            fft.getKSpace() = data;
            fft.invTransform();
            baseCoeff = fft.getRSpace();
        }
        else {
            initBaseCoeff();
            if (dataDim != data.getShape()) {
                dataDim = data.getShape();
                lu = DenseLU<T, true>(makeMatrixM());
                const VectorND<Tc> ones(lu.getOrder(), Tc(1));
                solverBuffer = lu.inv() * ones;
            }
            FFT<Tc, 3> fft(dataDim, PlanFlag::Estimate);
            auto lambdas = fft.getKSpace().flatten();
            lambdas = lu.inv() * data.flatten();

            const Tc average = lambdas.sum() / solverBuffer.sum();
            lambdas -= solverBuffer * average;
            fft.rawInvTransform();

            const Index3D baseDim = getBaseDim();
            baseCoeff.forND([this, &fft](Tc& coeff, Index3D index) {
                PeriodIndex3D pIndex(index, dataDim);
                pIndex.normalize();
                coeff *= fft.getRSpace()[pIndex];
            });
            baseCoeff[0, 0, 0] = average;
        }
    }

    template<Scalar T>
    T BrillouInterpolate<T>::interpolateFEM(Vec3D kPoint, const DenseTensor<T, 3>& data) const {
        using ElementType = CuboidLinear<Tr>;
        using IndexArray = ElementType::IndexArray;
        for (size_t i = 0; i < Dim; ++i) {
            kPoint[i] -= floor(kPoint[i]);
            assert(!kPoint[i].isNegative());
        }

        const Index3D shape = data.getShape();
        Vec3D globalPos{};
        for (int i = 0; i < 3; ++i)
            globalPos[i] = kPoint[i] * Tr(shape[i]);
        const Index3D index0{size_t(double(globalPos[0])), size_t(double(globalPos[1])), size_t(double(globalPos[2]))};
        const Vec3D point1{Tr(index0[0]), Tr(index0[1]), Tr(index0[2])};
        const Vec3D point2{Tr(index0[0] + 1), Tr(index0[1] + 1), Tr(index0[2] + 1)};
        const IndexArray nodeArr{
            ElementType::LeftFrontBottom,
            ElementType::LeftFrontTop,
            ElementType::LeftBehindBottom,
            ElementType::LeftBehindTop,
            ElementType::RightFrontBottom,
            ElementType::RightFrontTop,
            ElementType::RightBehindBottom,
            ElementType::RightBehindTop
        };
        const auto element = ElementType(point1, point2, nodeArr);
        const Vec3D localPos = element.toLocalPos(globalPos);

        T result = 0;
        size_t localNodeIndex = 0;
        for (int x = 0; x < 2; ++x) {
            for (int y = 0; y < 2; ++y) {
                for (int z = 0; z < 2; ++z) {
                    Index3D index{(index0[0] + x) % shape[0], (index0[1] + y) % shape[1], index0[2] + z};
                    const Tr factor = ElementType::baseFunc(nodeArr[localNodeIndex], localPos);
                    result += factor * data(index);
                    localNodeIndex += 1;
                }
            }
        }
        return result;
    }

    template<Scalar T>
    template<Matrix M>
    M BrillouInterpolate<T>::interpolateHermiteMatrix(Vec3D qPoint, const ArrayND<M, 3>& matrixGrid) {
        const size_t order = matrixGrid[0, 0, 0].getRow();
        M result(order, order, T(0));
        DenseTensor<Tc, 3> buffer(matrixGrid.getShape());
        for (size_t c = 0; c < order; ++c) {
            for (size_t r = 0; r < order; ++r) {
                buffer.forND([r, c, &matrixGrid](auto& elem, Index3D index) {
                    elem = matrixGrid(index)[r, c];
                });
                interpolate(buffer);
                const auto value = this->operator()(qPoint);
                result[r, c] = value;
            }
        }
        return result;
    }

    template<Scalar T>
    void BrillouInterpolate<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        baseCoeff.swap(obj.baseCoeff);
        lattice.swap(obj.lattice);
        lu.swap(obj.lu);
        solverBuffer.swap(solverBuffer);
        smoothFactor1.swap(obj.smoothFactor1);
        smoothFactor2.swap(obj.smoothFactor2);
        dataDim.swap(obj.dataDim);
    }

    template<Scalar T>
    void BrillouInterpolate<T>::initBaseCoeff() {
        const auto kernel = [this](Vec3D r, Index3D index) {
            const Tr r2 = r.squaredNorm();
            const bool isGammaPoint = r2 < std::numeric_limits<Tr>::epsilon();
            if (isGammaPoint)
                baseCoeff[index] = Tr(0);
            else
                baseCoeff[index] = reciprocal(Tr(1) + (smoothFactor1 + smoothFactor2 * r2) * r2);
        };
        PeriodicCell<Tr, 3>::template forLatticeCloud<true>(kernel, lattice, getBaseDim());
    }
    /**
     * Matrix M as defined in [1]
     */
    template<Scalar T>
    auto BrillouInterpolate<T>::makeMatrixM() const -> CoeffMatrixM {
        const size_t numDataPoint = dataDim[0] * dataDim[1] * dataDim[2];
        CoeffMatrixM matrixM(numDataPoint, numDataPoint);
        for (size_t r = 0; r < numDataPoint; ++r) {
            const PeriodIndex3D period_r(r, dataDim);
            for (size_t c = r; c < numDataPoint; ++c) {
                const PeriodIndex3D period_c(c, dataDim);
                const auto kIndex = Index3D(period_r + period_c);
                Tc elem(0);
                baseCoeff.forND([this, kIndex, &elem](const auto& coeff, Index3D index) {
                    Tr phase(0);
                    for (size_t i = 0; i < Dim; ++i)
                        phase += Tr(index[i] * kIndex[i]) / Tr(dataDim[i]);
                    phase *= Tr(2 * M_PI);
                    const auto factor = cos(phase);
                    elem += coeff * factor;
                });
                matrixM[r, c] = matrixM[c, r] = elem;
            }
        }
        return matrixM;
    }
}
