/*
 * Copyright 2024-2025 Weibo He.
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

namespace Physica {
    /**
     * Bohn-Huang(BH) projector impose Bohn-Huang and Huang invariance to force constant matrix as introduced in [1].
     * 
     * Reference:
     * [1] npj Comput. Mater. 8, 236 (2022); https://doi.org/10.1038/s41524-022-00920-6
     */
    template<Scalar T>
    class BHProjector final : public FCProjector<T> {
        using This = BHProjector<T>;
        using Base = FCProjector<T>;
    public:
        using typename Base::VectorType;
        using Base::Dim;
        using MDCellType = MDCell<T>;
        constexpr static unsigned int NumIndependentHuangCond = 15;
    private:
        using Base::superSize;
        using Base::numDOF;
        using Base::transBases;
        Array<VectorType> rotBases;
        Array<VectorType> huangBases;
    public:
        BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell);
        BHProjector(const This&) = default;
        BHProjector(This&&) noexcept = default;
        ~BHProjector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::projectSwap;
        T projectRot(VectorType& v) const;
        T projectHuang(VectorType& v) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumForceConsts;
    private:
        VectorType makeRotBase(size_t dof, unsigned int direction, const MDCellType& unitCell) const;
        VectorType makeHuangBase(unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta, const MDCellType& unitCell) const;
        void initHuangBases(const MDCellType& unitCell);
    };

    template<Scalar T>
    BHProjector<T>::BHProjector(Index3D superSize_, size_t numDOF_, const MDCellType& unitCell) {
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

    template<Scalar T>
    T BHProjector<T>::projectRot(VectorType& v) const {
        T maxAbsDot = 0;
        for (const auto& rotBase : rotBases) {
            const T dot = rotBase * v;
            maxAbsDot = std::max(maxAbsDot, abs(dot));
            v -= dot * rotBase;
        }
        return maxAbsDot;
    }

    template<Scalar T>
    T BHProjector<T>::projectHuang(VectorType& v) const {
        T maxAbsDot = 0;
        for (const auto& huangBase : huangBases) {
            const T dot = huangBase * v;
            maxAbsDot = std::max(maxAbsDot, abs(dot));
            v -= dot * huangBase;
        }
        return maxAbsDot;
    }

    template<Scalar T>
    void BHProjector<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
    }

    template<Scalar T>
    BHProjector<T>::VectorType BHProjector<T>::makeRotBase(
            size_t dof, unsigned int direction, const MDCellType& unitCell) const {
        const unsigned int dir1 = (direction + 1) % Dim;
        const unsigned int dir2 = (direction + 2) % Dim;
        VectorType result(getNumForceConsts(), 0);
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<T>::index1DTo5D(numDOF, superSize, i);
            if (index5D[3] == dof) {
                const size_t atom = index5D[4] / Dim;
                const unsigned int dir = index5D[4] % Dim;

                Vector3D<T> pos = unitCell.getPos().row(atom);
                for (unsigned int j = 0; j < Dim; ++j) {
                    const ssize_t index_j = index5D[j];
                    const ssize_t superSize_j = superSize[j];
                    const bool isOnWignerSeitzBoundary = (superSize_j % 2 == 0) && (index_j == superSize_j / 2);
                    if (!isOnWignerSeitzBoundary) {
                        const T factor = T(index_j > superSize_j / 2 ? index_j - superSize_j : index_j);
                        pos += unitCell.getLattice().row(j) * factor;
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

    template<Scalar T>
    BHProjector<T>::VectorType BHProjector<T>::makeHuangBase(
            unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta, const MDCellType& unitCell) const {
        assert(alpha < Dim && beta < Dim && gamma < Dim && delta < Dim && "[Error]: Invalid dimension");
        const auto& lattice = unitCell.getLattice();
        const auto& pos = unitCell.getPos();
        VectorType result(getNumForceConsts());
        for (size_t i = 0; i < result.getLength(); ++i) {
            const auto index5D = FCSwapVector<T>::index1DTo5D(numDOF, superSize, i);
            const unsigned int dir1 = index5D[3] % Dim;
            const unsigned int dir2 = index5D[4] % Dim;
            const Vector3D<T> pos1 = pos.row(index5D[3] / Dim);
            Vector3D<T> pos2 = pos.row(index5D[4] / Dim);
            for (unsigned int j = 0; j < Dim; ++j) {
                const ssize_t index_j = index5D[j];
                const ssize_t superSize_j = superSize[j];
                const bool isOnWignerSeitzBoundary = (superSize_j % 2 == 0) && (index_j == superSize_j / 2);
                if (!isOnWignerSeitzBoundary) {
                    const T factor = T(index_j > superSize_j / 2 ? index_j - superSize_j : index_j);
                    pos2 += lattice.row(j) * factor;
                }
            }

            T elem = 0;
            if (dir1 == alpha && dir2 == beta)
                elem += (pos2[gamma] - pos1[gamma]) * (pos2[delta] - pos1[delta]);
            if (dir1 == gamma && dir2 == delta)
                elem -= (pos2[alpha] - pos1[alpha]) * (pos2[beta] - pos1[beta]);
            result[i] = elem;
        }
        return result;
    }

    template<Scalar T>
    void BHProjector<T>::initHuangBases(const MDCellType& unitCell) {
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
