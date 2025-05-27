/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/PeriodIndex3D.h"
#include "Physica/Core/Math/Optimization/OptimizationImpl/QuadraticSearch.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "PhononSolver.h"

namespace Physica {
    /**
     * References:
     * [1] Dario Alfè PHON: A program to calculate phonons using the small displacement method [J]. Computer Physics Communications, 2009, 180(12), 2622-2633 (DOI: 10.1016/j.cpc.2009.03.010)
     */
    template<Scalar T>
    class FrozenPhonon : public PhononSolver<T> {
        using This = FrozenPhonon<T>;
        using Base = PhononSolver<T>;
    public:
        using typename Base::ComplexType;
        using typename Base::FFT3D;
        using typename Base::MDCellType;
        using typename Base::PositionMatrix;
        using typename Base::RSpaceFCGrid;
        constexpr static unsigned int Dim = Base::Dim;

        T displace;
    public:
        FrozenPhonon(MDCellType unitCell, Index3D superSize, T displace_);
        FrozenPhonon(const FrozenPhonon&) = default;
        FrozenPhonon(FrozenPhonon&&) noexcept = default;
        ~FrozenPhonon() = default;
        /* Operators */
        FrozenPhonon& operator=(FrozenPhonon obj) noexcept;
        /* Operations */
        template<class ForceModel>
        [[nodiscard]] RSpaceFCGrid makeForceConstants(ForceModel& model) const;
        template<class ForceModel>
        void optimize(const ForceModel& unitCellModel,
                      const ForceModel& superCellModel,
                      T minFreq,
                      unsigned int numInterpolateStep,
                      unsigned int maxNumStep);
        void swap(FrozenPhonon& __restrict obj) noexcept;
        /* Getters */
        using Base::getUnitCell;
        using Base::getNumUnitCellAtom;
        using Base::getUnitCellDOF;
        using Base::getSuperSize;
        using Base::getNumCell;
        /* Static members */
        static void read(RSpaceFCGrid& rSpaceFC, const H5Loc& loc, const char* name);
        static void write(const RSpaceFCGrid& rSpaceFC, H5Loc& loc, const char* name);
    private:
        PositionMatrix makeWignerSeitzRadius() const;
        ArrayND<DenseMatrix<T>, 3> makeWignerSeitzWeights() const;
        T calcWignerSeitzWeight(const Vector3D<T> r, const PositionMatrix& wignerSeitzRadius) const;
    };

    template<Scalar T>
    FrozenPhonon<T>::FrozenPhonon(MDCellType unitCell, Index3D superSize, T displace_)
            : Base(std::move(unitCell), superSize), displace(displace_) {}

    template<Scalar T>
    FrozenPhonon<T>& FrozenPhonon<T>::operator=(FrozenPhonon obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    template<class ForceModel>
    auto FrozenPhonon<T>::makeForceConstants(ForceModel& model) const -> RSpaceFCGrid {
        const size_t unitCellDOF = getUnitCellDOF();
        const T factor = -reciprocal(displace);
        const Index3D superSize = Base::getSuperSize();
        const MDCellType superCell = Base::getUnitCell().template makeSuperCell<ExtendCellOption::CellMajor>(superSize);

        RSpaceFCGrid result(getSuperSize(), unitCellDOF, unitCellDOF, T(0));
        PositionMatrix pos = superCell.getPos();
        for (size_t major = 0; major < unitCellDOF; ++major) {
            T& toDisplace = pos(major / Dim, major % Dim);
            const T copy = toDisplace;
            toDisplace += displace;
            VectorND<T> forceConst =
                    model.template force<Sequential>(MDCellType(superCell.getLattice(), pos, superCell.getMassVec())) * factor;
            toDisplace = copy;

            for (size_t cell = 0; cell < getNumCell(); ++cell) {
                const Index3D index = PeriodIndex3D(cell, superSize);
                for (size_t minor = 0; minor < unitCellDOF; ++minor)
                    result(index).refFromMajorMinor(major, minor) = forceConst[cell * unitCellDOF + minor];
            }
        }
        return result;
    }

    template<Scalar T>
    template<class ForceModel>
    void FrozenPhonon<T>::optimize(
            const ForceModel& unitCellModel,
            const ForceModel& superCellModel,
            T minFreq,
            unsigned int numInterpolateStep,
            unsigned int maxNumStep) {
        assert(minFreq.isPositive() && "[Error]: Min frequency must be positive");
        for (unsigned int step = 0; step < maxNumStep; ++step) {
            const auto fcMatrixes = makeForceConstants(superCellModel);
            auto fcMatrix = fcMatrixes(0, 0, 0);
            Base::toDynamicMatrix(fcMatrix);
            const auto eigen = Base::diagonalize(fcMatrix);
            const T freq = Base::makeFreq(eigen)[0];
            const bool isImagFreqSmallEnough = -freq < minFreq;
            if (isImagFreqSmallEnough)
                return;

            const auto eigenvectors = Base::makeEigenVectors(eigen);
            const T e0 = unitCellModel.potentialV(Base::getUnitCell());
            const T e1 = unitCellModel.potentialV(Base::shiftAtom(eigenvectors.col(0), displace));
            const T e2 = unitCellModel.potentialV(Base::shiftAtom(eigenvectors.col(0), -displace));
            QuadraticSearch<T> optimizer({0, displace, -displace}, {e0, e1, e2});
            for (unsigned int i = 0; i < numInterpolateStep; ++i) {
                optimizer.step([this, &unitCellModel, &eigenvectors](T dist) {
                    return unitCellModel.potentialV(Base::shiftAtom(eigenvectors.col(0), dist));
                });
            }
            Base::setUnitCell(Base::shiftAtom(eigenvectors.col(0), optimizer.getOptimizedX()));
        }
    }

    template<Scalar T>
    void FrozenPhonon<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        displace.swap(obj.displace);
    }

    template<Scalar T>
    FrozenPhonon<T>::PositionMatrix
    FrozenPhonon<T>::makeWignerSeitzRadius() const {
        constexpr int OneSideDim = 2; // 2 is enough to generate a Wigner-Seitz cell
        constexpr int TwoSideDim = 2 * OneSideDim + 1;
        constexpr size_t ResultSize = TwoSideDim * TwoSideDim * TwoSideDim - 1;
        const Index3D superSize = Base::getSuperSize();
        PositionMatrix wignerSeitzRadius(ResultSize, 3);
        size_t index = 0;
        for (int i = -OneSideDim; i <= OneSideDim; ++i) {
            for (int j = -OneSideDim; j <= OneSideDim; ++j) {
                for (int k = -OneSideDim; k <= OneSideDim; ++k) {
                    const Vector3D<T> factor{T(i * superSize[0]), T(j * superSize[1]), T(k * superSize[2])};
                    auto row = wignerSeitzRadius.row(index);
                    row = Base::getUnitCell().getLattice().transpose() * factor;
                    const bool isNotGammaPoint = row.squaredNorm() > std::numeric_limits<T>::epsilon();
                    index += isNotGammaPoint;
                }
            }
        }
        assert(index == ResultSize && "[Error]: Wrong result");
        return wignerSeitzRadius;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    void FrozenPhonon<T>::read(RSpaceFCGrid& rSpaceFC, const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        unsigned char superSize[Dim];
        auto attr = group.openAttribute("SuperSize");
        attr.read(H5::PredType::NATIVE_UCHAR, superSize);
        const auto gridDim = Index3D{superSize[0], superSize[1], superSize[2]};
        rSpaceFC.resize(gridDim);

        rSpaceFC.forND([&group](auto& fc, Index3D index) {
            fc.read(group, std::format("{}_{}_{}", index[0], index[1], index[2]).c_str());
        });
    }

    template<Scalar T>
    void FrozenPhonon<T>::write(const RSpaceFCGrid& rSpaceFC, H5Loc& loc, const char* name) {
        auto group = loc.openGroup(name);
        /* Write superSize */ {
            const auto shape = rSpaceFC.getShape();
            unsigned char superSize[Dim];
            for (unsigned int i = 0; i < Dim; ++i) {
                assert(shape[i] <= std::numeric_limits<unsigned char>::max() && "[Error]: Unexpected large super cell");
                superSize[i] = shape[i];
            }
            group.writeAttr("SuperSize", superSize);
        }

        rSpaceFC.forND([&group, &rSpaceFC](const auto& fc, Index3D index) {
            fc.write(group, std::format("{}_{}_{}", index[0], index[1], index[2]).c_str());
        });
    }
#endif
    template<Scalar T>
    ArrayND<DenseMatrix<T>, 3> FrozenPhonon<T>::makeWignerSeitzWeights() const {
        const auto wignerSeitzRadius = makeWignerSeitzRadius();
        const Index3D superSize = Base::getSuperSize();
        const Index3D gridDim{4 * superSize[0] + 1, 4 * superSize[1] + 1, 4 * superSize[2] + 1};
        const size_t numAtom = getNumUnitCellAtom();
        ArrayND<DenseMatrix<T>, 3> result(gridDim, numAtom, numAtom);
        forND(gridDim, [this, superSize, numAtom, &result, &wignerSeitzRadius](Index3D index) {
            const auto& unitCell = Base::getUnitCell();
            const Vector3D<T> factor{T(index[0]) - T(2 * superSize[0]),
                                  T(index[1]) - T(2 * superSize[1]),
                                  T(index[2]) - T(2 * superSize[2])};
            const Vector3D<T> r0 = unitCell.getLattice().transpose() * factor;
            auto& mat = result(index);
            for (size_t c = 0; c < result.getCol(); ++c) {
                for (size_t r = 0; r < result.getRow(); ++r) {
                    const Vector3D<T> r1 = r0 + unitCell.getPos().row(r) - unitCell.getPos().row(c);
                    mat(r, c) = calcWignerSeitzWeight(r1, wignerSeitzRadius);
                }
            }
        });
        return result;
    }

    template<Scalar T>
    T FrozenPhonon<T>::calcWignerSeitzWeight(
            const Vector3D<T> r, const PositionMatrix& wignerSeitzRadius) const {
        constexpr double precision = 1E-5;
        size_t count = 1;
        for (size_t i = 0; i < wignerSeitzRadius.getRow(); ++i) {
            const auto radius = wignerSeitzRadius.row(i);
            const T dot = r * radius;
            const T halfSquaredNorm = radius.squaredNorm() * T(0.5);
            const bool isOnBoundary = scalarNear(dot, halfSquaredNorm, precision);
            if (isOnBoundary) [[unlikely]]
                count += 1;
            else {
                const bool isOutsideWignerSeitzCell = dot > halfSquaredNorm;
                if (isOutsideWignerSeitzCell)
                    return T(0);
            }
        }
        return reciprocal(T(count));
    }
}
