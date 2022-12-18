/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    /**
     * References:
     * [1] M. Matsumoto, T. Yagasaki, H. Tanaka. GenIce: Hydrogen-Disordered Ice Generator.[J]. J. Comput. Chem. 2017, DOI: 10.1002/jcc.25077
     */
    class RandIce {
        using ScalarType = typename CrystalCell::ScalarType;
        using PositionMatrix = typename CrystalCell::PositionMatrix;
        using Vector3D = Vector<ScalarType, 3>;
        constexpr static double BondLengthOH = PhyConst<AU>::angstormToBohr(1);

        CrystalCell initialCell;
        ScalarType maxDistOO;
        ScalarType maxDistOH;
        Utils::Array<bool> isHydrogenOccupied;
        Utils::Array<unsigned char> numHydrogenRequired;
    public:
        RandIce(CrystalCell initialCell_, ScalarType maxDistOO_, ScalarType maxDistOH_);
        RandIce(const RandIce&) = default;
        RandIce(RandIce&&) noexcept = default;
        ~RandIce() = default;
        /* Operators */
        RandIce& operator=(RandIce obj) noexcept;
        /* Operations */
        template<class RandomGenerator> CrystalCell run(RandomGenerator& gen);
        void swap(RandIce& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return initialCell.getNumParticle() / 3U; }
    private:
        void prepareRun();
        void searchDanglingH(PositionMatrix& pos);
        Utils::Array<size_t> findOInRadius(size_t indexO, ScalarType radius) const;
        Utils::Array<size_t> findBondedHInRadius(size_t indexO, ScalarType radius) const;
        template<class RandomGenerator> size_t makeRandEmptyO(RandomGenerator& gen) const;
        template<class RandomGenerator> size_t makeRandFreeH(size_t indexO, RandomGenerator& gen) const;
        void fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH);
        size_t countFreeH(const Utils::Array<size_t>& hIndexes) const;
        void searchForPairs(PositionMatrix& pos);
        size_t getIndexToPair() const;
        template<class RandomGenerator> void randUninitializedH(PositionMatrix& pos, RandomGenerator& gen);
        /* Getters */
        [[nodiscard]] size_t getEndIndexH() const noexcept { return getNumMolecule() * 2U; }
        [[nodiscard]] size_t getStartIndexO() const noexcept { return getEndIndexH(); }
        [[nodiscard]] size_t getEndIndexO() const noexcept { return initialCell.getNumParticle(); }
        [[nodiscard]] bool isFinished() const noexcept;
        /* Static members */
        template<class RandomGenerator>
        static Vector<ScalarType, 3> randUnitVector(RandomGenerator& gen);
    };

    template<class RandomGenerator>
    CrystalCell RandIce::run(RandomGenerator& gen) {
        PositionMatrix pos = initialCell.getPos();
        prepareRun();
        searchDanglingH(pos);
        while (!isFinished()) {
            const size_t randO = makeRandEmptyO(gen);
            const size_t randH = makeRandFreeH(randO, gen);
            fetchHydrogen(pos, randO, randH);
            searchForPairs(pos);
        }
        randUninitializedH(pos, gen);
        CrystalCell result(initialCell.getLattice(), std::move(pos), initialCell.getAtomicNumbers(), CrystalCell::Type::Cartesian);
        result.normalizeCartesianCell();
        return result;
    }

    template<class RandomGenerator>
    size_t RandIce::makeRandEmptyO(RandomGenerator& gen) const {
        std::uniform_int_distribution<size_t> dist(0, getNumMolecule() - 1);
        size_t result = dist(gen);
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto hInRange = findBondedHInRadius(result, maxDistOO);
            const size_t numFreeH = countFreeH(hInRange);
            const bool isOxygenEmpty = numHydrogenRequired[result] != 0;
            if (isOxygenEmpty && numFreeH != 0)
                break;
            result = (result + 1U) % getNumMolecule();
        }
        assert(numHydrogenRequired[result] != 0);
        return result;
    }

    template<class RandomGenerator>
    size_t RandIce::makeRandFreeH(size_t indexO, RandomGenerator& gen) const {
        assert(indexO < getNumMolecule());
        const auto hInRange = findBondedHInRadius(indexO, maxDistOO);
        size_t randLogicIndex;
        /* Rand index */ {
            const size_t numFreeH = countFreeH(hInRange);
            assert(numFreeH > numHydrogenRequired[indexO]);
            std::uniform_int_distribution<size_t> dist(0, numFreeH - 1);
            randLogicIndex = dist(gen);
        }
        size_t logicIndex = 0;
        for (size_t physicalIndex = 0; physicalIndex < hInRange.getLength(); ++physicalIndex) {
            const bool isFree = isHydrogenOccupied[hInRange[physicalIndex]] == false;
            if (isFree && logicIndex == randLogicIndex)
                return hInRange[physicalIndex];
            logicIndex += isFree;
        }
        [[maybe_unused]] constexpr bool ShouldNotReachHere = false;
        assert(ShouldNotReachHere);
        return hInRange[0];
    }

    template<class RandomGenerator>
    void RandIce::randUninitializedH(PositionMatrix& pos, RandomGenerator& gen) {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            while (numHydrogenRequired[i] != 0) {
                auto row = pos.row(i * 2U + (2U - numHydrogenRequired[i]));
                row = initialCell.getPos().row(i + getStartIndexO()).asVector() + randUnitVector(gen) * ScalarType(BondLengthOH);
                numHydrogenRequired[i] -= 1;
            }
        }
    }

    template<class RandomGenerator>
    typename RandIce::Vector3D RandIce::randUnitVector(RandomGenerator& gen) {
        std::uniform_real_distribution dist{};
        const ScalarType theta(dist(gen) * M_PI);
        const ScalarType phi(dist(gen) * M_PI * 2);
        Vector<ScalarType, 3> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
        return result;
    }
}
