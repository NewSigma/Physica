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
        constexpr static unsigned int NumIndependentHuangCond = 15;
    private:
        using Base::superSize;
        using Base::numDOF;
        using Base::transBases;
        Utils::Array<VectorType> rotBases;
        Utils::Array<VectorType> huangBases;
    public:
        BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell);
        BHProjector(const BHProjector&) = default;
        BHProjector(BHProjector&&) noexcept = default;
        ~BHProjector() = default;
        /* Operators */
        BHProjector& operator=(BHProjector obj) noexcept { swap(obj); return obj; }
        /* Operations */
        using Base::projectSwap;
        ScalarType projectRot(VectorType& v) const;
        ScalarType projectHuang(VectorType& v) const;
        void swap(BHProjector& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumForceConsts;
    private:
        VectorType makeRotBase(size_t dof, unsigned int direction, const MDCellType& unitCell) const;
        VectorType makeHuangBase(unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta, const MDCellType& unitCell) const;
        void initHuangBases(const MDCellType& unitCell);
    };

    template<class ScalarType>
    BHProjector<ScalarType>::BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell) {
        superSize = std::move(superSize_);
        numDOF = numDOF_;
        const size_t numBase = numDOF * Dim;
        transBases.reserve(numBase);
        rotBases.reserve(numBase);
        for (size_t i = 0; i < numBase; ++i) {
            VectorType transBase = Base::makeTransBase(i / Dim, i % Dim);
            transBase.toUnit();
            transBases.grow(std::move(transBase));

            VectorType rotBase = makeRotBase(i / Dim, i % Dim, unitCell);
            rotBase.toUnit();
            rotBases.grow(std::move(rotBase));
        }
        initHuangBases(unitCell);
    }

    template<class ScalarType>
    ScalarType BHProjector<ScalarType>::projectRot(VectorType& v) const {
        ScalarType maxAbsDot = 0;
        for (const auto& rotBase : rotBases) {
            const ScalarType dot = rotBase * v;
            maxAbsDot = std::max(maxAbsDot, abs(dot));
            v -= dot * rotBase;
        }
        return maxAbsDot;
    }

    template<class ScalarType>
    ScalarType BHProjector<ScalarType>::projectHuang(VectorType& v) const {
        ScalarType maxAbsDot = 0;
        for (const auto& huangBase : huangBases) {
            const ScalarType dot = huangBase * v;
            maxAbsDot = std::max(maxAbsDot, abs(dot));
            v -= dot * huangBase;
        }
        return maxAbsDot;
    }

    template<class ScalarType>
    void BHProjector<ScalarType>::swap(BHProjector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
    }

    template<class ScalarType>
    typename BHProjector<ScalarType>::VectorType BHProjector<ScalarType>::makeRotBase(
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

    template<class ScalarType>
    typename BHProjector<ScalarType>::VectorType BHProjector<ScalarType>::makeHuangBase(
            unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta, const MDCellType& unitCell) const {
        assert(alpha < Dim && beta < Dim && gamma < Dim && delta < Dim && "[Error]: Invalid dimension");
        const auto& lattice = unitCell.getLattice();
        const auto& pos = unitCell.getPos();
        VectorType result(getNumForceConsts());
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<ScalarType>::index1DTo5D(numDOF, superSize, i);
            const unsigned int dir1 = index5D[3] % Dim;
            const unsigned int dir2 = index5D[4] % Dim;
            const Vector<ScalarType, 3> pos1 = pos.row(index5D[3] / Dim);
            Vector<ScalarType, 3> pos2 = pos.row(index5D[4] / Dim);
            for (unsigned int j = 0; j < Dim; ++j) {
                const ssize_t index_j = index5D[j];
                const ssize_t superSize_j = superSize[j];
                const bool isOnWignerSeitzBoundary = (superSize_j % 2 == 0) && (index_j == superSize_j / 2);
                if (!isOnWignerSeitzBoundary) {
                    const ScalarType factor = ScalarType(index_j > superSize_j / 2 ? index_j - superSize_j : index_j);
                    pos2 += lattice.row(j).asVector() * factor;
                }
            }

            ScalarType elem = 0;
            if (dir1 == alpha && dir2 == beta)
                elem += (pos2[gamma] - pos1[gamma]) * (pos2[delta] - pos1[delta]);
            if (dir1 == gamma && dir2 == delta)
                elem -= (pos2[alpha] - pos1[alpha]) * (pos2[beta] - pos1[beta]);
            result[i] = elem;
        }
        return result;
    }

    template<class ScalarType>
    void BHProjector<ScalarType>::initHuangBases(const MDCellType& unitCell) {
        huangBases.reserve(NumIndependentHuangCond);
        huangBases.grow(makeHuangBase(0, 0, 0, 1, unitCell));
        huangBases.grow(makeHuangBase(0, 0, 0, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 0, 1, 1, unitCell));
        huangBases.grow(makeHuangBase(0, 0, 1, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 0, 2, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 1, 0, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 1, 1, 1, unitCell));
        huangBases.grow(makeHuangBase(0, 1, 1, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 1, 2, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 2, 1, 1, unitCell));
        huangBases.grow(makeHuangBase(0, 2, 1, 2, unitCell));
        huangBases.grow(makeHuangBase(0, 2, 2, 2, unitCell));
        huangBases.grow(makeHuangBase(1, 1, 1, 2, unitCell));
        huangBases.grow(makeHuangBase(1, 1, 2, 2, unitCell));
        huangBases.grow(makeHuangBase(1, 2, 2, 2, unitCell));
        for (auto& base : huangBases)
            base.toUnit();
    }
}
