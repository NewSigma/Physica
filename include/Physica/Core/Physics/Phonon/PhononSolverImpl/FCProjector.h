/*
 * Copyright 2024 WeiBo He.
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

#include "FCSwapVector.h"

namespace Physica::Core {
    /**
     * Force constant(FC) projector utilizes the projection method to incorporate symmetry invariance as introduced in [1].
     * 
     * References:
     * [1] Mounet, N. Structural, vibrational and thermodynamic properties of carbon allotropes from ﬁrst-principles: diamond, graphite, and nanotubes. Master’s thesis, Massachusetts Institute of Technology (2005); http://hdl.handle.net/1721.1/33400
     */
    template<class ScalarType>
    class FCProjector {
    public:
        constexpr static unsigned int Dim = 3;
        using VectorType = Vector<ScalarType>;
        using RSpaceFCMat = DenseMatrix<ScalarType>;
        using RSpaceFCGrid = GridStorage<RSpaceFCMat>;
        using Index3D = typename RSpaceFCGrid::Index3D;
    protected:
        Index3D superSize;
        size_t numDOF;
        Utils::Array<VectorType> transBases;
    public:
        FCProjector() = default;
        FCProjector(Index3D superSize_, size_t numDOF_);
        FCProjector(const FCProjector&) = default;
        FCProjector(FCProjector&&) noexcept = default;
        ~FCProjector() = default;
        /* Operators */
        FCProjector& operator=(FCProjector obj) noexcept { swap(obj); return obj; }
        /* Operations */
        ScalarType projectSwap(VectorType& v) const;
        ScalarType projectTrans(VectorType& v) const;
        VectorType toVector(const RSpaceFCGrid& fcGrid) const;
        void toGrid(const VectorType& fcVector, RSpaceFCGrid& fcGrid) const;
        void swap(FCProjector& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Index3D getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSize[0] * superSize[1] * superSize[2]; }
        [[nodiscard]] size_t getNumDOF() const noexcept { return numDOF; }
        [[nodiscard]] size_t getNumUnitCellAtom() const noexcept { return numDOF / Dim; }
        [[nodiscard]] size_t getNumSuperCellAtom() const noexcept { return getNumUnitCellAtom() * getNumCell(); }
        [[nodiscard]] inline size_t getNumForceConsts() const noexcept;
    protected:
        VectorType makeTransBase(size_t dof, unsigned int direction) const;
    };

    template<class ScalarType>
    FCProjector<ScalarType>::FCProjector(Index3D superSize_, size_t numDOF_)
            : superSize(std::move(superSize_)), numDOF(numDOF_) {
        const size_t numBase = numDOF * Dim;
        transBases.reserve(numBase);
        for (size_t i = 0; i < numBase; ++i) {
            VectorType transBase = makeTransBase(i / Dim, i % Dim);
            transBase.toUnit();
            transBases.grow(std::move(transBase));
        }
    }

    template<class ScalarType>
    ScalarType FCProjector<ScalarType>::projectSwap(VectorType& v) const {
        ScalarType maxAbsDot = 0;
        GridBase::forIndexInGrid(superSize, [this, &v, &maxAbsDot](Index3D index) {
            ScalarType maxAbsDotInner = 0;
            for (size_t i = 0; i < numDOF; ++i) {
                for (size_t j = i; j < numDOF; ++j) {
                    const bool isSelf = i == j && index[0] == 0 && index[1] == 0 && index[2] == 0;
                    if (isSelf) [[unlikely]]
                        continue;

                    const FCSwapVector<ScalarType> swapBase(getSuperSize(), index, numDOF, i, j);
                    const ScalarType dot = swapBase * v;
                    maxAbsDotInner = std::max(maxAbsDotInner, abs(dot));
                    v -= dot * swapBase;
                }
            }
            maxAbsDot = std::max(maxAbsDot, maxAbsDotInner);
        });
        return maxAbsDot;
    }

    template<class ScalarType>
    ScalarType FCProjector<ScalarType>::projectTrans(VectorType& v) const {
        ScalarType maxAbsDot = 0;
        for (const auto& transBase : transBases) {
            const ScalarType dot = transBase * v;
            maxAbsDot = std::max(maxAbsDot, abs(dot));
            v -= dot * transBase;
        }
        return maxAbsDot;
    }

    template<class ScalarType>
    typename FCProjector<ScalarType>::VectorType FCProjector<ScalarType>::toVector(const RSpaceFCGrid& fcGrid) const {
        assert(fcGrid.getDim() == superSize && "[Error]: Cell sizes do not match");
        assert(fcGrid(0, 0, 0).getRow() == numDOF && "[Error]: DOFs do not match");
        VectorType result(getNumForceConsts());
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            const Index3D cellIndex{index5D[0], index5D[1], index5D[2]};
            result[i] = fcGrid(cellIndex).calcFromMajorMinor(index5D[3], index5D[4]);
        }
        return result;
    }

    template<class ScalarType>
    void FCProjector<ScalarType>::toGrid(const VectorType& fcVector, RSpaceFCGrid& fcGrid) const {
        assert(fcVector.getLength() == getNumForceConsts() && "[Error]: This is not a force constants vector");
        assert(fcGrid.getDim() == superSize && "[Error]: Cell sizes do not match");
        assert(fcGrid(0, 0, 0).getRow() == numDOF && "[Error]: DOFs do not match");
        for (size_t i = 0; i < fcVector.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            const Index3D cellIndex{index5D[0], index5D[1], index5D[2]};
            fcGrid(cellIndex).refFromMajorMinor(index5D[3], index5D[4]) = fcVector[i];
        }
    }

    template<class ScalarType>
    void FCProjector<ScalarType>::swap(FCProjector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        superSize.swap(obj.superSize);
        std::swap(numDOF, obj.numDOF);
        transBases.swap(obj.transBases);
    }

    template<class ScalarType>
    inline size_t FCProjector<ScalarType>::getNumForceConsts() const noexcept {
        return numDOF * numDOF * superSize[0] * superSize[1] * superSize[2];
    }

    template<class ScalarType>
    typename FCProjector<ScalarType>::VectorType FCProjector<ScalarType>::makeTransBase(size_t dof, unsigned int direction) const {
        VectorType result(getNumForceConsts(), 0);
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            if ((index5D[3] == dof) && (index5D[4] % Dim == direction))
                result[i] = ScalarType(1);
        }
        return result;
    }
}
