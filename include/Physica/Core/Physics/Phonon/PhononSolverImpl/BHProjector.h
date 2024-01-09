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

#include "FCProjector.h"

namespace Physica::Core {
    /**
     * Bohn-Huang(BH) projector impose Bohn-Huang and Huang invariance to force constant matrix as introduced in [1].
     * 
     * Reference:
     * [1] npj Comput. Mater. 8, 236 (2022); https://doi.org/10.1038/s41524-022-00920-6
     */
    template<class ScalarType>
    class BHProjector final : public FCProjector<ScalarType> {
        using Base = FCProjector<ScalarType>;
    public:
        using typename Base::VectorType;
        using typename Base::Index3D;
        using Base::Dim;
        using MDCellType = MDCell<ScalarType>;
    private:
        using Base::superSize;
        using Base::numDOF;
        using Base::transBases;
        Utils::Array<VectorType> rotBases;
    public:
        BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell);
        BHProjector(const BHProjector&) = default;
        BHProjector(BHProjector&&) noexcept = default;
        ~BHProjector() = default;
        /* Operators */
        BHProjector& operator=(BHProjector obj) noexcept { swap(obj); return obj; }
        /* Operations */
        using Base::projectSwap;
        void projectRot(VectorType& v) const;
        void swap(BHProjector& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumForceConsts;
    private:
        VectorType makeFCRotVector(size_t dof, unsigned int direction, const MDCellType& unitCell) const;
    };

    template<class ScalarType>
    BHProjector<ScalarType>::BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell) {
        assert(pos.getRow() == getNumSuperCellAtom() && "[Error]: Invalid pos");
        superSize = std::move(superSize_);
        numDOF = numDOF_;
        const size_t numBase = numDOF * Dim;
        transBases.reserve(numBase);
        rotBases.reserve(numBase);
        for (size_t i = 0; i < numBase; ++i) {
            VectorType transBase = Base::makeFCTransVector(i / Dim, i % Dim);
            projectSwap(transBase);
            transBases.grow(std::move(transBase));

            VectorType rotBase = makeFCRotVector(i / Dim, i % Dim, unitCell);
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
    void BHProjector<ScalarType>::projectRot(VectorType& v) const {
        for (const auto& rotBase : rotBases)
            v -= (rotBase * v) * rotBase;
    }

    template<class ScalarType>
    void BHProjector<ScalarType>::swap(BHProjector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
    }

    template<class ScalarType>
    typename BHProjector<ScalarType>::VectorType BHProjector<ScalarType>::makeFCRotVector(
            size_t dof, unsigned int direction, const MDCellType& unitCell) const {
        const unsigned int dir1 = (direction + 1) % Dim;
        const unsigned int dir2 = (direction + 2) % Dim;
        VectorType result(getNumForceConsts(), 0);
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            if (index5D[3] == dof) {
                const size_t atom = index5D[4] / Dim;
                const unsigned int dir = index5D[4] % Dim;

                Vector<ScalarType, 3> pos = unitCell.getPos().row(atom);
                for (unsigned int j = 0; j < Dim; ++j) {
                    const ssize_t index_j = index5D[j];
                    const ssize_t superSize_j = superSize[j];
                    const bool isOnWignerSeitzBoundary = (superSize_j % 2 == 0) && (index_j == superSize_j / 2);
                    if (!isOnWignerSeitzBoundary) {
                        const ScalarType factor = ScalarType(index_j > superSize_j / 2 ? index_j - superSize_j : index_j);
                        pos += unitCell.getLattice().row(j).asVector() * factor;
                    }
                }
                if (dir == dir1)
                    result[i] = -pos[dir2];
                else if (dir == dir2)
                    result[i] = pos[dir1];
            }
        }
        return result;
    }
}
