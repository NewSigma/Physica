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
#include "Physica/Core/Physics/RandIce.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    RandIce::RandIce(CrystalCell initialCell_, ScalarType maxDistOO_, ScalarType maxDistOH_)
            : initialCell(std::move(initialCell_))
            , maxDistOO(maxDistOO_)
            , maxDistOH(maxDistOH_) {
        assert(initialCell.getNumParticle() % 3U == 0);
        isHydrogenOccupied.resize(getEndIndexH());
        numHydrogenRequired.resize(getNumMolecule());
        if (initialCell.getType() == CrystalCell::Type::Direct)
            initialCell.toCartesian();
    }

    RandIce& RandIce::operator=(RandIce obj) noexcept {
        swap(obj);
        return *this;
    }

    void RandIce::swap(RandIce& obj) noexcept {
        initialCell.swap(obj.initialCell);
        std::swap(maxDistOO, obj.maxDistOO);
        std::swap(maxDistOH, obj.maxDistOH);
        isHydrogenOccupied.swap(obj.isHydrogenOccupied);
        numHydrogenRequired.swap(obj.numHydrogenRequired);
    }

    void RandIce::prepareRun() {
        for (auto& elem : isHydrogenOccupied)
            elem = false;
        for (auto& elem : numHydrogenRequired)
            elem = 2U;
    }

    void RandIce::searchDanglingH(PositionMatrix& pos) {
        const size_t startIndexO = getStartIndexO();
        const size_t endIndexO = getEndIndexO();
        for (size_t i = startIndexO; i < endIndexO; ++i) {
            const auto hInRadius = findBondedHInRadius(i - startIndexO, maxDistOO);
            unsigned char numDanglingH = 0;
            for (size_t j = 0; j < getEndIndexH(); ++j) {
                const bool isInRange = initialCell.minDistVector(i, j).squaredNorm() < square(maxDistOH);
                const bool noHydrogenBond = std::find(hInRadius.begin(), hInRadius.end(), j) == hInRadius.end();
                const bool isDanglingH = isInRange && noHydrogenBond;
                if (isDanglingH) {
                    numDanglingH += 1;
                    isHydrogenOccupied[j] = true;
                    auto r = pos.row((i - startIndexO) * 2U + numDanglingH);
                    r = initialCell.getPos().row(j).asVector();
                }
            }
            numHydrogenRequired[i - startIndexO] -= numDanglingH;
            if (numDanglingH > 2U)
                throw std::invalid_argument("[Error]: Bad initial structure");
        }
    }

    Utils::Array<size_t> RandIce::findOInRadius(size_t indexO, ScalarType radius) const {
        Utils::Array<size_t> result{};
        const ScalarType squaredRadiusO = square(radius);
        for (size_t i = getStartIndexO(); i < getEndIndexO(); ++i) {
            const size_t shiftIndexO = getStartIndexO() + indexO;
            if (i == shiftIndexO)
                continue;
            const bool isInRadius = initialCell.minDistVector(shiftIndexO, i).squaredNorm() < squaredRadiusO;
            if (isInRadius)
                result.append(i);
        }
        return result;
    }

    Utils::Array<size_t> RandIce::findBondedHInRadius(size_t indexO, ScalarType radius) const {
        const auto otherO = findOInRadius(indexO, radius);
        const ScalarType squaredRadiusH = square(radius * 0.5);
        Utils::Array<size_t> result{};
        for (size_t oxygen : otherO) {
            const size_t shiftIndexO = getStartIndexO() + indexO;
            const Vector3D delta = initialCell.minDistVector(shiftIndexO, oxygen);
            const Vector3D middle = initialCell.getPos().row(shiftIndexO).asVector() + delta * ScalarType(0.5);
            ScalarType minSquaredDist = std::numeric_limits<ScalarType>::max();
            size_t indexH = 0;
            for (size_t i = 0; i < getEndIndexH(); ++i) {
                const ScalarType squaredDist = initialCell.minDistVector(middle, i).squaredNorm();
                if (squaredDist < minSquaredDist) {
                    minSquaredDist = squaredDist;
                    indexH = i;
                }
            }
            result.append(indexH);
        }
        return result;
    }

    std::pair<size_t, size_t> RandIce::findHydrogenInMolecule(size_t indexO) const {
        ScalarType minSquaredDist1 = std::numeric_limits<ScalarType>::max();
        ScalarType minSquaredDist2 = std::numeric_limits<ScalarType>::max();
        size_t indexH1 = 0;
        size_t indexH2 = 1;
        for (size_t i = 0; i < getEndIndexH(); ++i) {
            const ScalarType squaredDist = initialCell.minDistVector(i, getStartIndexO() + indexO).squaredNorm();
            if (squaredDist < minSquaredDist1) {
                indexH1 = i;
                minSquaredDist1 = squaredDist;
            }
            else if (squaredDist < minSquaredDist2) {
                indexH2 = i;
                minSquaredDist2 = squaredDist;
            }

            if (minSquaredDist1 < minSquaredDist2) {
                std::swap(minSquaredDist1, minSquaredDist2);
                std::swap(indexH1, indexH2);
            }
        }
        return {indexH1, indexH2};
    }

    void RandIce::fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH) {
        assert(!isHydrogenOccupied[indexH]);
        assert(numHydrogenRequired[indexO] > 0);
        
        const size_t shiftIndexO = getStartIndexO() + indexO;
        Vector3D delta = initialCell.minDistVector(shiftIndexO, indexH);
        delta.toUnit();
        delta *= ScalarType(BondLengthOH);
        auto row = pos.row(indexO * 2U + (2U - numHydrogenRequired[indexO]));
        row = initialCell.getPos().row(shiftIndexO).asVector() + delta;

        isHydrogenOccupied[indexH] = true;
        numHydrogenRequired[indexO] -= 1;
    }

    size_t RandIce::countFreeH(const Utils::Array<size_t>& hIndexes) const {
        size_t numFreeH = 0;
        for (auto h : hIndexes)
            numFreeH += isHydrogenOccupied[h] == false;
        return numFreeH;
    }

    void RandIce::searchForPairs(PositionMatrix& pos) {
        size_t indexO = getIndexToPair();
        bool isValidIndex = indexO != getNumMolecule();
        while (isValidIndex) {
            const auto hInRange = findBondedHInRadius(indexO, maxDistOO);
            for (auto h : hInRange) {
                if (!isHydrogenOccupied[h])
                    fetchHydrogen(pos, indexO, h);
            }
            indexO = getIndexToPair();
            isValidIndex = indexO != getNumMolecule();
        }
    }

    size_t RandIce::getIndexToPair() const {
        for (size_t indexO = 0; indexO < getNumMolecule(); ++indexO) {
            const auto hInRange = findBondedHInRadius(indexO, maxDistOO);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH > 0 && numFreeH == numHydrogenRequired[indexO])
                return indexO;
        }
        return getNumMolecule();
    }

    bool RandIce::isFinished() const noexcept {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto hInRange = findBondedHInRadius(i, maxDistOO);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH != 0 && numHydrogenRequired[i] != 0)
                return false;
        }
        return true;
    }
}
