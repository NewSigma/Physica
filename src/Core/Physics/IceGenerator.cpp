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
#include <algorithm>
#include "Physica/Core/Physics/IceGenerator.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    IceGenerator::IceGenerator(CrystalCell initialCell_, ScalarType maxDistOO_, ScalarType maxDistOH_)
            : initialCell(std::move(initialCell_))
            , maxDistOO(maxDistOO_)
            , maxDistOH(maxDistOH_) {
        assert(initialCell.getNumParticle() % 3U == 0);
        isHydrogenOccupied.resize(getEndIndexH());
        numHydrogenRequired.resize(getNumMolecule());
        if (initialCell.getType() == CrystalCell::Type::Direct)
            initialCell.toCartesian();
    }

    IceGenerator& IceGenerator::operator=(IceGenerator obj) noexcept {
        swap(obj);
        return *this;
    }

    Utils::Array<CrystalCell> IceGenerator::exhaust() {
        Utils::Array<CrystalCell> result{};
        prepareRun();
        searchDanglingH();
        exhaustImpl(0, initialCell.getPos(), result);
        return result;
    }

    void IceGenerator::swap(IceGenerator& obj) noexcept {
        initialCell.swap(obj.initialCell);
        std::swap(maxDistOO, obj.maxDistOO);
        std::swap(maxDistOH, obj.maxDistOH);
        isHydrogenOccupied.swap(obj.isHydrogenOccupied);
        numHydrogenRequired.swap(obj.numHydrogenRequired);
    }

    void IceGenerator::prepareRun() {
        for (auto& elem : isHydrogenOccupied)
            elem = false;
        for (auto& elem : numHydrogenRequired)
            elem = 2U;
    }

    void IceGenerator::searchDanglingH() {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto pair = findHydrogenInMolecule(i);
            const auto hInRadius = findBondedH(i);
            const bool isDanglingH1 = std::find(hInRadius.begin(), hInRadius.end(), pair.first) == hInRadius.end();
            const bool isDanglingH2 = std::find(hInRadius.begin(), hInRadius.end(), pair.second) == hInRadius.end();
            isHydrogenOccupied[pair.first] = isDanglingH1;
            isHydrogenOccupied[pair.second] = isDanglingH2;
            numHydrogenRequired[i] -= static_cast<unsigned char>(isDanglingH1) + static_cast<unsigned char>(isDanglingH2);
        }
    }

    Utils::Array<size_t> IceGenerator::findOInRadius(size_t indexO, ScalarType radius) const {
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

    Utils::Array<size_t> IceGenerator::findBondedH(size_t indexO) const {
        const auto range = CrystalCell::estimateRange(initialCell.getLattice(), maxDistOO);
        Utils::Array<size_t> result{};
        CrystalCell::forCellInRange(range, initialCell.getLattice(), [this, indexO, &result](Vector3D delta) {
            const ScalarType squaredRadiusH = square(maxDistOO * 0.5);
            const ScalarType squaredRadiusO = square(maxDistOO);
            const size_t shiftIndexO = getStartIndexO() + indexO;

            for (size_t i = getStartIndexO(); i < getEndIndexO(); ++i) {
                const Vector3D otherO = initialCell.getPos().row(i) + delta;
                const ScalarType r2 = (initialCell.getPos().row(shiftIndexO) - otherO).squaredNorm();
                const bool isNotSelf = std::numeric_limits<ScalarType>::epsilon() < r2;
                const bool isInRange = r2 < squaredRadiusO;
                if (isNotSelf && isInRange) {
                    const Vector3D middle = (otherO + initialCell.getPos().row(shiftIndexO)) * ScalarType(0.5);
                    for (size_t j = 0; j < getEndIndexH(); ++j) {
                        const bool isHInRange = initialCell.minDistVector(middle, j).squaredNorm() < squaredRadiusH;
                        const bool isOccupied = std::find(result.cbegin(), result.cend(), j) != result.cend();
                        if (isHInRange && !isOccupied)
                            result.append(j);
                    }
                }
            }
        });
        return result;
    }

    Utils::Array<size_t> IceGenerator::findFreeBondedHInRadius(size_t indexO) const {
        Utils::Array<size_t> result{};
        const auto hInRange = findBondedH(indexO);
        for (auto h : hInRange)
            if (!isHydrogenOccupied[h])
                result.append(h);
        return result;
    }

    std::pair<size_t, size_t> IceGenerator::findHydrogenInMolecule(size_t indexO) const {
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

    void IceGenerator::fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH) {
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

    size_t IceGenerator::countFreeH(const Utils::Array<size_t>& hIndexes) const {
        size_t numFreeH = 0;
        for (auto h : hIndexes)
            numFreeH += isHydrogenOccupied[h] == false;
        return numFreeH;
    }

    void IceGenerator::searchForPairs(PositionMatrix& pos) {
        size_t indexO = getIndexToPair();
        bool isValidIndex = indexO != getNumMolecule();
        while (isValidIndex) {
            const auto hInRange = findBondedH(indexO);
            for (auto h : hInRange) {
                if (!isHydrogenOccupied[h])
                    fetchHydrogen(pos, indexO, h);
            }
            indexO = getIndexToPair();
            isValidIndex = indexO != getNumMolecule();
        }
    }

    size_t IceGenerator::getIndexToPair() const {
        for (size_t indexO = 0; indexO < getNumMolecule(); ++indexO) {
            const auto hInRange = findBondedH(indexO);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH > 0 && numFreeH == numHydrogenRequired[indexO])
                return indexO;
        }
        return getNumMolecule();
    }

    void IceGenerator::exhaustImpl(size_t stackDepth, const PositionMatrix& pos, Utils::Array<CrystalCell>& result) {
        const bool recursionStop = stackDepth == getNumMolecule();
        if (recursionStop) {
            CrystalCell cell({initialCell.getLattice(), pos, CrystalCell::Type::Cartesian}, initialCell.getAtomicNumbers());
            cell.normalize();
            result.append(std::move(cell));
            return;
        }

        const unsigned char numNeedH = numHydrogenRequired[stackDepth];
        if (numNeedH == 0) {
            exhaustImpl(stackDepth + 1, pos, result);
            return;
        }

        const auto freeH = findFreeBondedHInRadius(stackDepth);
        if (freeH.getLength() < numNeedH)
            return;

        for (size_t j = 0; j < freeH.getLength() - 1; ++j) {
            PositionMatrix copy = pos;
            fetchHydrogen(copy, stackDepth, freeH[j]);
            if (numNeedH == 2U) {
                for (size_t k = j + 1; k < freeH.getLength(); ++k) {
                    PositionMatrix copy1 = copy;
                    fetchHydrogen(copy1, stackDepth, freeH[k]);
                    exhaustImpl(stackDepth + 1, copy1, result);

                    numHydrogenRequired[stackDepth] += 1;
                    isHydrogenOccupied[freeH[k]] = false;
                }
            }
            else
                exhaustImpl(stackDepth + 1, copy, result);

            numHydrogenRequired[stackDepth] += 1;
            isHydrogenOccupied[freeH[j]] = false;
        }
    }

    bool IceGenerator::isFinished() const noexcept {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto hInRange = findBondedH(i);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH != 0 && numHydrogenRequired[i] != 0)
                return false;
        }
        return true;
    }
}
