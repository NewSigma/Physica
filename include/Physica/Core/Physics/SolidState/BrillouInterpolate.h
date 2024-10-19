/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations/LUSolver.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/CuboidLinear.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica::Core {
    /**
     * References:
     * [1] R.N. Euwema, D.J. Stukel, T.C. Collins, J.S. DeWitt, D.G. Shankland, Phys. Rev. 178 (1969) 1419–1423. http://dx.doi.org/10.1103/PhysRev.178.1419.
     */
    template<class ScalarType>
    class BrillouInterpolate {
        constexpr static unsigned int Dim = 3;
        constexpr static bool isComplex = ScalarType::isComplex;
        using Index3D = typename GridBase::Index3D;
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
        using LatticeMatrix = typename PeriodicCell<RealType, Dim>::LatticeMatrix;
        using CoeffMatrixM = DenseMatrix<ComplexType, MatrixOption::Row | MatrixOption::Vector>;
        using Vector3D = Vector<RealType, 3>;

        RSpaceGrid<ComplexType> baseCoeff;
        LatticeMatrix lattice;
        LUSolver<ComplexType, CoeffMatrixM::Option> solver;
        Vector<ComplexType> solverBuffer;
        RealType smoothFactor1;
        RealType smoothFactor2;
        Index3D dataDim;
    public:
        BrillouInterpolate(Index3D baseDim, LatticeMatrix unitcell, RealType smoothFactor1_, RealType smoothFactor2_);
        BrillouInterpolate(const BrillouInterpolate& obj) = default;
        BrillouInterpolate(BrillouInterpolate&& obj) noexcept = default;
        ~BrillouInterpolate() = default;
        /* Operators */
        BrillouInterpolate& operator=(BrillouInterpolate obj) noexcept;
        [[nodiscard]] ScalarType operator()(Vector3D kPoint);
        /* Operations */
        RealType calcRoughness() const;
        void interpolate(const RSpaceGrid<ScalarType>& data);
        ScalarType interpolateFEM(Vector3D kPoint, const RSpaceGrid<ScalarType>& data) const;
        template<class MatrixType>
        MatrixType interpolateHermiteMatrix(Vector3D qPoint, const GridStorage<MatrixType>& matrixGrid);
        void swap(BrillouInterpolate& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Index3D getBaseDim() const noexcept { return baseCoeff.getDim(); }
    private:
        void initBaseCoeff();
        CoeffMatrixM makeMatrixM() const;
    };

    template<class ScalarType>
    BrillouInterpolate<ScalarType>::BrillouInterpolate(Index3D baseDim, LatticeMatrix unitcell, RealType smoothFactor1_, RealType smoothFactor2_)
            : baseCoeff(baseDim)
            , lattice(unitcell)
            , smoothFactor1(smoothFactor1_)
            , smoothFactor2(smoothFactor2_)
            , dataDim({0, 0, 0}) {
        RealType minSquaredNorm = std::numeric_limits<RealType>::max();
        for (int i = 0; i < 3; ++i)
            minSquaredNorm = std::min(minSquaredNorm, lattice.row(i).squaredNorm());
        smoothFactor1 /= minSquaredNorm;
        smoothFactor2 /= square(minSquaredNorm);
    }

    template<class ScalarType>
    ScalarType BrillouInterpolate<ScalarType>::operator()(Vector3D kPoint) {
        const Index3D baseDim = getBaseDim();
        ScalarType result(0);
        GridBase::forIndexInGrid(baseDim, [this, kPoint, baseDim, &result](Index3D index) {
            RealType phase(0);
            for (size_t dim = 0; dim < Dim; ++dim)
                phase += RealType(index[dim]) * kPoint[dim];
            phase *= RealType(2 * M_PI);
            const auto factor = cos(phase);

            if constexpr (!isComplex) {
                PeriodIndex3D pIndex(index, baseDim);
                if (pIndex.isInReducedK())
                    result += (baseCoeff(index) * factor).real();
                else
                    result += (baseCoeff(pIndex.toReducedK()).conjugate() * factor).real();
            }
            else {
                (void)baseDim;
                result += baseCoeff(index) * factor;
            }
        });
        return result;
    }

    template<class ScalarType>
    typename BrillouInterpolate<ScalarType>::RealType BrillouInterpolate<ScalarType>::calcRoughness() const {
        RealType result = 0;
        const auto kernel = [this, &result](Vector3D r, Index3D index) {
            const RealType r2 = r.squaredNorm();
            const bool isGammaPoint = r2 < std::numeric_limits<RealType>::epsilon();
            if (isGammaPoint)
                return;
            const RealType r4 = square(r2);
            result += baseCoeff(index).squaredNorm() * (RealType(1) + smoothFactor1 * r2 + smoothFactor2 * r4);
        };
        GridBase::forPointIndexInGrid<true, decltype(kernel)>(getBaseDim(), lattice, kernel);
        return result;
    }

    template<class ScalarType>
    void BrillouInterpolate<ScalarType>::interpolate(const RSpaceGrid<ScalarType>& data) {
        assert(dataDim[0] <= getBaseDim()[0] && "[Error]: Not enough base function, interpolation is overdetermined");
        assert(dataDim[1] <= getBaseDim()[1] && "[Error]: Not enough base function, interpolation is overdetermined");
        assert(dataDim[2] <= getBaseDim()[2] && "[Error]: Not enough base function, interpolation is overdetermined");
        const bool useNormalFFT = dataDim[0] == getBaseDim()[0] && dataDim[1] == getBaseDim()[1] && dataDim[2] == getBaseDim()[2];
        if (useNormalFFT) {
            dataDim = data.getDim();
            FFT<ComplexType, 3> fft(dataDim, PlanFlag::Estimate);
            fft.getKSpace() = data;
            fft.invTransform();
            baseCoeff = fft.getRSpace();
        }
        else {
            initBaseCoeff();
            if (dataDim != data.getDim()) {
                dataDim = data.getDim();
                solver.decomposition(makeMatrixM());
                const Vector<ComplexType> ones(solver.getOrder(), ComplexType(1));
                solverBuffer = solver.solve(ones);
            }
            const auto v1 = solver.solve(data.flatten());
            const ComplexType average = v1.sum() / solverBuffer.sum();

            FFT<ComplexType, 3> fft(dataDim, PlanFlag::Estimate);
            auto lambdas = fft.getKSpace().flatten();
            lambdas = v1 - solverBuffer * average;
            fft.rawInvTransform();

            const Index3D baseDim = getBaseDim();
            GridBase::forIndexInGrid(baseDim, [this, &fft](Index3D rIndex) {
                PeriodIndex3D pIndex(rIndex, dataDim);
                pIndex.normalize();
                baseCoeff(rIndex) *= fft.getRSpace()(pIndex);
            });
            baseCoeff(0, 0, 0) = average;
        }
    }

    template<class ScalarType>
    ScalarType BrillouInterpolate<ScalarType>::interpolateFEM(Vector3D kPoint, const RSpaceGrid<ScalarType>& data) const {
        using ElementType = CuboidLinear<RealType>;
        using IndexArray = typename ElementType::IndexArray;
        for (size_t i = 0; i < Dim; ++i) {
            kPoint[i] -= floor(kPoint[i]);
            assert(!kPoint[i].isNegative());
        }

        const Index3D dataDim = data.getDim();
        Vector3D globalPos{};
        for (int i = 0; i < 3; ++i)
            globalPos[i] = kPoint[i] * RealType(dataDim[i]);
        const Index3D index0{size_t(double(globalPos[0])), size_t(double(globalPos[1])), size_t(double(globalPos[2]))};
        const Vector3D point1{RealType(index0[0]), RealType(index0[1]), RealType(index0[2])};
        const Vector3D point2{RealType(index0[0] + 1), RealType(index0[1] + 1), RealType(index0[2] + 1)};
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
        const Vector3D localPos = element.toLocalPos(globalPos);

        ScalarType result = 0;
        size_t localNodeIndex = 0;
        for (int x = 0; x < 2; ++x) {
            for (int y = 0; y < 2; ++y) {
                for (int z = 0; z < 2; ++z) {
                    Index3D index{(index0[0] + x) % dataDim[0], (index0[1] + y) % dataDim[1], index0[2] + z};
                    const RealType factor = ElementType::baseFunc(nodeArr[localNodeIndex], localPos);
                    result += factor * data(index);
                    localNodeIndex += 1;
                }
            }
        }
        return result;
    }

    template<class ScalarType>
    template<class MatrixType>
    MatrixType BrillouInterpolate<ScalarType>::interpolateHermiteMatrix(Vector3D qPoint, const GridStorage<MatrixType>& matrixGrid) {
        const Index3D gridDim = matrixGrid.getDim();
        const size_t order = matrixGrid(0, 0, 0).getRow();
        MatrixType result(order, order, ScalarType(0));
        RSpaceGrid<ComplexType> buffer(gridDim);
        for (size_t c = 0; c < order; ++c) {
            for (size_t r = 0; r < order; ++r) {
                GridBase::forIndexInGrid(gridDim, [r, c, &buffer, &matrixGrid](Index3D index) {
                    buffer(index) = matrixGrid(index)(r, c);
                });
                interpolate(buffer);
                const auto value = this->operator()(qPoint);
                result(r, c) = value;
            }
        }
        return result;
    }

    template<class ScalarType>
    void BrillouInterpolate<ScalarType>::swap(BrillouInterpolate& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        baseCoeff.swap(obj.baseCoeff);
        lattice.swap(obj.lattice);
        solver.swap(obj.solver);
        solverBuffer.swap(solverBuffer);
        smoothFactor1.swap(obj.smoothFactor1);
        smoothFactor2.swap(obj.smoothFactor2);
        dataDim.swap(obj.dataDim);
    }

    template<class ScalarType>
    void BrillouInterpolate<ScalarType>::initBaseCoeff() {
        const auto kernel = [this](Vector3D r, Index3D index) {
            const RealType r2 = r.squaredNorm();
            const bool isGammaPoint = r2 < std::numeric_limits<RealType>::epsilon();
            if (isGammaPoint)
                baseCoeff(index) = RealType(0);
            else
                baseCoeff(index) = reciprocal(RealType(1) + (smoothFactor1 + smoothFactor2 * r2) * r2);
        };
        GridBase::forPointIndexInGrid<RealType, true, decltype(kernel)>(getBaseDim(), lattice, kernel);
    }
    /**
     * Matrix M as defined in [1]
     */
    template<class ScalarType>
    typename BrillouInterpolate<ScalarType>::CoeffMatrixM BrillouInterpolate<ScalarType>::makeMatrixM() const {
        const auto baseDim = getBaseDim();
        const size_t numDataPoint = dataDim[0] * dataDim[1] * dataDim[2];
        CoeffMatrixM matrixM(numDataPoint, numDataPoint);
        for (size_t r = 0; r < numDataPoint; ++r) {
            const PeriodIndex3D period_r(r, dataDim);
            for (size_t c = r; c < numDataPoint; ++c) {
                const PeriodIndex3D period_c(c, dataDim);
                const auto kIndex = Index3D(period_r + period_c);
                ComplexType elem(0);
                GridBase::forIndexInGrid(baseDim, [this, kIndex, &elem](Index3D rIndex) {
                    RealType phase(0);
                    for (size_t i = 0; i < Dim; ++i)
                        phase += RealType(rIndex[i] * kIndex[i]) / RealType(dataDim[i]);
                    phase *= RealType(2 * M_PI);
                    const auto factor = cos(phase);
                    elem += baseCoeff(rIndex) * factor;
                });
                matrixM(r, c) = matrixM(c, r) = elem;
            }
        }
        return matrixM;
    }
}
