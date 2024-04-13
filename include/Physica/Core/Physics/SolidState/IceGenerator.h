/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Utils/TestHelper.h"

namespace Physica::Core {
    /**
     * References:
     * [1] M. Matsumoto, T. Yagasaki, H. Tanaka. GenIce: Hydrogen-Disordered Ice Generator.[J]. J. Comput. Chem. 2017, DOI: 10.1002/jcc.25077
     */
    template<class ScalarType>
    class IceGenerator {
        using PositionMatrix = typename PeriodicCell<ScalarType, 3>::PositionMatrix;
        using Vector3D = Vector<ScalarType, 3>;
        using CrystalCellType = CrystalCell<ScalarType>;
        constexpr static double BondLengthOH = PhyConst<AU>::angstormToBohr(1);

        CrystalCellType initialCell;
        ScalarType maxDistOO;
        ScalarType maxDistOH;
        Utils::Array<bool> isHydrogenOccupied;
        Utils::Array<unsigned char> numHydrogenRequired;
    public:
        IceGenerator(ScalarType maxDistOO_, ScalarType maxDistOH_);
        IceGenerator(CrystalCellType initialCell_, ScalarType maxDistOO_, ScalarType maxDistOH_);
        IceGenerator(const IceGenerator&) = default;
        IceGenerator(IceGenerator&&) noexcept = default;
        ~IceGenerator() = default;
        /* Operators */
        IceGenerator& operator=(IceGenerator obj) noexcept;
        /* Operations */
        Utils::Array<CrystalCellType> exhaust();
        template<class RandomGenerator> CrystalCellType makeRand(RandomGenerator& gen);
        template<class RandomGenerator> CrystalCellType makeDefects(unsigned int numDefect, RandomGenerator& gen) const;
        template<class RandomGenerator> Utils::Array<size_t> randRing(RandomGenerator& gen) const;
        CrystalCellType makeRingMove(const Utils::Array<size_t>& ring, PositionMatrix& momentumMat) const;
        void swap(IceGenerator& __restrict obj) noexcept;
        /* Setters */
        void setInitialCell(CrystalCellType cell);
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return initialCell.getNumParticle() / 3U; }
    private:
        PositionMatrix prepareRun();
        void searchDanglingH(PositionMatrix& pos);
        Utils::Array<size_t> findOInRadius(size_t indexO, ScalarType radius) const;
        Utils::Array<size_t> findBondedH(size_t indexMolecule) const;
        Utils::Array<size_t> findFreeBondedHInRadius(size_t indexO) const;
        std::pair<size_t, size_t> findHydrogenInMolecule(size_t indexMolecule) const;
        size_t findHydrogenBetweenO(size_t indexO1, size_t indexO2) const;
        template<class RandomGenerator> size_t makeRandEmptyO(RandomGenerator& gen) const;
        template<class RandomGenerator> size_t makeRandFreeH(size_t indexO, RandomGenerator& gen) const;
        void fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH);
        size_t countFreeH(const Utils::Array<size_t>& hIndexes) const;
        void searchForPairs(PositionMatrix& pos);
        size_t getIndexToPair() const;
        template<class RandomGenerator> void randUninitializedH(PositionMatrix& pos, RandomGenerator& gen);
        void exhaustImpl(size_t stackDepth, const PositionMatrix& pos, Utils::Array<CrystalCellType>& result);
        /* Getters */
        [[nodiscard]] size_t getEndIndexH() const noexcept { return getNumMolecule() * 2U; }
        [[nodiscard]] size_t getStartIndexO() const noexcept { return getEndIndexH(); }
        [[nodiscard]] size_t getEndIndexO() const noexcept { return initialCell.getNumParticle(); }
        [[nodiscard]] bool isFinished() const noexcept;
        /* Static members */
        static Vector3D rotate(ScalarType angle, Vector3D axis, Vector3D target);
        template<class RandomGenerator>
        static Vector<ScalarType, 3> randUnitVector(RandomGenerator& gen);

        friend class ::Physica::Test;
    };
}

#include "IceGeneratorImpl.h"
