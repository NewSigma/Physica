/*
 * Copyright 2022-2024 Weibo He.
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

#include <algorithm>
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica {
    /**
     * References:
     * [1] M. Matsumoto, T. Yagasaki, H. Tanaka. GenIce: Hydrogen-Disordered Ice Generator.[J]. J. Comput. Chem. 2017, DOI: 10.1002/jcc.25077
     */
    template<Scalar T>
    class IceGenerator {
        using This = IceGenerator<T>;
        using PositionMatrix = PeriodicCell<T, 3>::PositionMatrix;
        using CrystalCellType = CrystalCell<T>;
        constexpr static double BondLengthOH = PhyConst<AU>::angstormToBohr(1);

        CrystalCellType initialCell;
        T maxDistOO;
        T maxDistOH;
        Array<bool> isHydrogenOccupied;
        Array<unsigned char> numHydrogenRequired;
    public:
        IceGenerator(T maxDistOO_, T maxDistOH_);
        IceGenerator(CrystalCellType initialCell_, T maxDistOO_, T maxDistOH_);
        IceGenerator(const This&) = default;
        IceGenerator(This&&) noexcept = default;
        ~IceGenerator() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        [[nodiscard]] Array<CrystalCellType> exhaust();
        template<RNG R>
        [[nodiscard]] CrystalCellType makeRand();
        template<RNG R>
        [[nodiscard]] CrystalCellType makeDefects(unsigned int numDefect) const;
        template<RNG R>
        [[nodiscard]] Array<size_t> randRing() const;
        [[nodiscard]] CrystalCellType makeRingMove(const Array<size_t>& ring, PositionMatrix& momentumMat) const;
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setInitialCell(CrystalCellType cell);
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return initialCell.getNumParticle() / 3U; }
    private:
        [[nodiscard]] PositionMatrix prepareRun();
        void searchDanglingH(PositionMatrix& pos);
        [[nodiscard]] Array<size_t> findOInRadius(size_t indexO, T radius) const;
        [[nodiscard]] Array<size_t> findBondedH(size_t indexMolecule) const;
        [[nodiscard]] Array<size_t> findFreeBondedHInRadius(size_t indexO) const;
        [[nodiscard]] std::pair<size_t, size_t> findHydrogenInMolecule(size_t indexMolecule) const;
        [[nodiscard]] size_t findHydrogenBetweenO(size_t indexO1, size_t indexO2) const;
        template<RNG R>
        [[nodiscard]] size_t makeRandEmptyO() const;
        template<RNG R>
        [[nodiscard]] size_t makeRandFreeH(size_t indexO) const;
        void fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH);
        [[nodiscard]] size_t countFreeH(const Array<size_t>& hIndexes) const;
        void searchForPairs(PositionMatrix& pos);
        [[nodiscard]] size_t getIndexToPair() const;
        template<RNG R> void randUninitializedH(PositionMatrix& pos);
        void exhaustImpl(size_t stackDepth, const PositionMatrix& pos, Array<CrystalCellType>& result);
        /* Getters */
        [[nodiscard]] size_t getEndIndexH() const noexcept { return getNumMolecule() * 2U; }
        [[nodiscard]] size_t getStartIndexO() const noexcept { return getEndIndexH(); }
        [[nodiscard]] size_t getEndIndexO() const noexcept { return initialCell.getNumParticle(); }
        [[nodiscard]] bool isFinished() const noexcept;
        /* Static members */
        static Vector3D<T> rotate(T angle, Vector3D<T> axis, Vector3D<T> target);
        template<RNG R>
        static Vector3D<T> randUnitVector();

        friend class ::Physica::Test;
    };
}

#include "IceGeneratorImpl.h"
