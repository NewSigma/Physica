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
        using MDCellType = MDCell<ScalarType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using Index3D = typename RSpaceFCGrid::Index3D;
    private:
        Index3D superSize;
        size_t numDOF;
        Utils::Array<VectorType> transBases;
        Utils::Array<VectorType> rotBases;
    public:
        FCProjector(Index3D superSize_, size_t numDOF_);
        FCProjector(Index3D superSize_, size_t numDOF_, const PositionMatrix& pos);
        FCProjector(const FCProjector&) = default;
        FCProjector(FCProjector&&) noexcept = default;
        ~FCProjector() = default;
        /* Operators */
        FCProjector& operator=(FCProjector obj) noexcept { swap(obj); return obj; }
        /* Operations */
        void projectSwap(VectorType& v) const;
        void projectTrans(VectorType& v) const;
        void projectRot(VectorType& v) const;
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
    private:
        VectorType makeFCTransVector(size_t dof, unsigned int direction) const;
        VectorType makeFCRotVector(size_t dof, unsigned int direction, const PositionMatrix& superPos) const;
    };

    template<class ScalarType>
    FCProjector<ScalarType>::FCProjector(Index3D superSize_, size_t numDOF_)
            : superSize(std::move(superSize_)), numDOF(numDOF_) {
        const size_t numBase = numDOF * Dim;
        transBases.reserve(numBase);
        for (size_t i = 0; i < numBase; ++i) {
            VectorType transBase = makeFCTransVector(i / Dim, i % Dim);
            projectSwap(transBase);
            transBases.grow(std::move(transBase));
        }

        for (size_t i = 0; i < numBase; ++i) {
            auto& base_i = transBases[i];
            assert(!base_i.norm().isZero() && "[Error]: Unexpected base degenerate");
            base_i.toUnit();

            for (size_t j = i + 1; j < numBase; ++j) {
                auto& base_j = transBases[j];
                base_j -= (base_i * base_j) * base_i;
            }
        }
    }

    template<class ScalarType>
    FCProjector<ScalarType>::FCProjector(Index3D superSize_, size_t numDOF_, const PositionMatrix& pos)
            : superSize(std::move(superSize_)), numDOF(numDOF_) {
        assert(pos.getRow() == getNumUnitCellAtom() && "[Error]: Invalid pos");
        const size_t numBase = numDOF * Dim;
        transBases.reserve(numBase);
        rotBases.reserve(numBase);
        for (size_t i = 0; i < numBase; ++i) {
            VectorType transBase = makeFCTransVector(i / Dim, i % Dim);
            projectSwap(transBase);
            transBases.grow(std::move(transBase));

            VectorType rotBase = makeFCRotVector(i / Dim, i % Dim, pos);
            rotBase.toUnit();
            projectSwap(rotBase);
            rotBases.grow(std::move(rotBase));
        }

        for (size_t i = 0; i < transBases.getLength(); ++i) {
            auto& base_i = transBases[i];
            assert(!base_i.norm().isZero() && "[Error]: Unexpected base degenerate");
            base_i.toUnit();

            for (size_t j = i + 1; j < transBases.getLength(); ++j) {
                auto& base_j = transBases[j];
                base_j -= (base_i * base_j) * base_i;
            }

            for (auto& rotBase : rotBases)
                rotBase -= (rotBase * base_i) * base_i;
        }

        for (size_t i = 0; i < rotBases.getLength(); ++i) {
            auto& base_i = rotBases[i];
            assert(!base_i.norm().isZero() && "[Error]: Unexpected base degenerate");
            base_i.toUnit();

            for (size_t j = i + 1; j < rotBases.getLength(); ++j) {
                auto& base_j = rotBases[j];
                base_j -= (base_i * base_j) * base_i;
            }
        }
    }

    template<class ScalarType>
    void FCProjector<ScalarType>::projectSwap(VectorType& v) const {
        GridBase::forIndexInGrid(superSize, [this, &v](Index3D index) {
            for (size_t i = 0; i < numDOF; ++i) {
                for (size_t j = i; j < numDOF; ++j) {
                    const bool isSelf = i == j && index[0] == 0 && index[1] == 0 && index[2] == 0;
                    if (isSelf) [[unlikely]]
                        continue;

                    const FCSwapVector<ScalarType> swapBase(getSuperSize(), index, numDOF, i, j);
                    v -= (swapBase * v) * swapBase;
                }
            }
        });
    }

    template<class ScalarType>
    void FCProjector<ScalarType>::projectTrans(VectorType& v) const {
        for (const auto& transBase : transBases)
            v -= (transBase * v) * transBase;
    }

    template<class ScalarType>
    void FCProjector<ScalarType>::projectRot(VectorType& v) const {
        for (const auto& rotBase : rotBases)
            v -= (rotBase * v) * rotBase;
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
    }

    template<class ScalarType>
    inline size_t FCProjector<ScalarType>::getNumForceConsts() const noexcept {
        return numDOF * numDOF * superSize[0] * superSize[1] * superSize[2];
    }

    template<class ScalarType>
    typename FCProjector<ScalarType>::VectorType FCProjector<ScalarType>::makeFCTransVector(size_t dof, unsigned int direction) const {
        VectorType result(getNumForceConsts(), 0);
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            if ((index5D[3] == dof) && (index5D[4] % Dim == direction))
                result[i] = ScalarType(1);
        }
        return result;
    }

    template<class ScalarType>
    typename FCProjector<ScalarType>::VectorType FCProjector<ScalarType>::makeFCRotVector(
            size_t dof, unsigned int direction, const PositionMatrix& pos) const {
        const unsigned int dir1 = (direction + 1) % Dim;
        const unsigned int dir2 = (direction + 2) % Dim;
        VectorType result(getNumForceConsts(), 0);
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            if (index5D[3] == dof) {
                const size_t atom = index5D[4] / Dim;
                const unsigned int dir = index5D[4] % Dim;
                if (dir == dir1)
                    result[i] = -pos(atom, dir2);
                else if (dir == dir2)
                    result[i] = pos(atom, dir1);
            }
        }
        return result;
    }
}
