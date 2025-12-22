/*
 * Copyright 2023-2024 Weibo He.
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
#include "IceGenerator.h"

namespace Physica {
    template<Scalar T>
    IceGenerator<T>::IceGenerator(T maxDistOO_, T maxDistOH_)
            : maxDistOO(maxDistOO_)
            , maxDistOH(maxDistOH_) {}

    template<Scalar T>
    IceGenerator<T>::IceGenerator(CrystalCellType initialCell_, T maxDistOO_, T maxDistOH_)
            : IceGenerator(maxDistOO_, maxDistOH_) {
        assert(initialCell.getNumParticle() % 3U == 0);
        setInitialCell(std::move(initialCell_));
    }

    template<Scalar T>
    auto IceGenerator<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    auto IceGenerator<T>::exhaust() -> Array<CrystalCellType> {
        Array<CrystalCellType> result{};
        PositionMatrix pos = prepareRun();
        searchDanglingH(pos);
        exhaustImpl(0, pos, result);
        return result;
    }

    template<Scalar T>
    template<RNG R>
    auto IceGenerator<T>::makeRand() -> CrystalCellType{
        PositionMatrix pos = prepareRun();
        searchDanglingH(pos);
        while (!isFinished()) {
            const size_t randO = makeRandEmptyO<R>();
            const size_t randH = makeRandFreeH<R>(randO);
            fetchHydrogen(pos, randO, randH);
            searchForPairs(pos);
        }
        randUninitializedH<R>(pos);
        CrystalCellType result(initialCell.getLattice(), std::move(pos), initialCell.getAtomicNumbers(), CrystalCellType::Type::Cartesian);
        result.normalize();
        return result;
    }

    template<Scalar T>
    template<RNG R>
    auto IceGenerator<T>::makeDefects(unsigned int numDefect) const -> CrystalCellType {
        assert(numDefect < getNumMolecule());
        PositionMatrix pos = initialCell.getPos();

        Array<size_t> permutation(getNumMolecule());
        for (size_t i = 0; i < permutation.getLength(); ++i)
            permutation[i] = i;

        auto& gen = R::getInstance();
        for (unsigned int i = 0; i < numDefect; ++i) {
            const size_t indexDefectO = [this, i]() {
                std::uniform_int_distribution<size_t> dist(0, getNumMolecule() - 1 - i);
                const size_t randIndex = dist(gen);
                size_t result = permutation[randIndex];
                std::swap(permutation[randIndex], permutation[getNumMolecule() - 1 - i]);
                return result;
            }();
            auto otherO = findOInRadius(indexDefectO, maxDistOO);
            const auto hydrogenInMolecular = findHydrogenInMolecule(indexDefectO);
            {
                std::uniform_int_distribution<size_t> dist(0, otherO.getLength() - 1);
                const size_t randIndexO = otherO[dist(gen)];
                const size_t randIndexH = (gen() % 2U == 0) ? hydrogenInMolecular.first : hydrogenInMolecular.second;
                auto row = pos.row(randIndexH);
                const Vector3D<T> delta = initialCell.minDistVector(getStartIndexO() + indexDefectO, randIndexO);
                row = pos.row(getStartIndexO() + indexDefectO) + delta * (T(BondLengthOH) / delta.norm());
            }
        }
        CrystalCellType result(initialCell.getLattice(), std::move(pos), initialCell.getAtomicNumbers(), CrystalCellType::Type::Cartesian);
        result.normalize();
        return result;
    }
    /**
     * Reference:
     * [1] J. Chem. Phys. 118, 9291 (2003); https://doi.org/10.1063/1.1568337
     */
    template<Scalar T>
    template<RNG R>
    Array<size_t> IceGenerator<T>::randRing() const {
        auto& gen = R::getInstance();
        Array<size_t> ring{};
        size_t ringStart = 0;
        {
            std::uniform_int_distribution<size_t> dist(0, getNumMolecule() - 1);
            const size_t randMolecule = dist(gen);
            ring.append(randMolecule);
        }

        bool isRingGenerated = false;
        while(!isRingGenerated) {
            const size_t lastMolecule = ring[ring.getLength() - 1];
            const auto pair = findHydrogenInMolecule(lastMolecule);
            const auto bondH = findBondedH(lastMolecule);
            const bool isBondedH1 = std::find(bondH.begin(), bondH.end(), pair.first) != bondH.end();
            const bool isBondedH2 = std::find(bondH.begin(), bondH.end(), pair.second) != bondH.end();
            size_t indexH = 0;
            if (isBondedH1 && isBondedH2)
                indexH = gen() % 2 == 0 ? pair.first : pair.second;
            else if (isBondedH1 || isBondedH2)
                indexH = isBondedH1 ? pair.first : pair.second;
            else
                return Array<size_t>{};
            
            size_t nextIndexO = getNumMolecule();
            /* Find next oxygen */ {
                const auto range = CrystalCellType::estimateRange(initialCell.getLattice(), maxDistOO);
                const T squaredRadiusO = square(maxDistOO);
                const size_t indexO = getStartIndexO() + lastMolecule;
                T minSquaredDist = std::numeric_limits<T>::max();
                for (size_t i = getStartIndexO(); i < getEndIndexO(); ++i) {
                    const T r2 = initialCell.minDistVector(i, indexO).squaredNorm();
                    const bool isNotSelf = r2 > T(std::numeric_limits<T>::epsilon());
                    const bool isInRange = r2 < squaredRadiusO;
                    if (isNotSelf && isInRange) {
                        const T squaredDist = initialCell.minDistVector(indexH, i).squaredNorm();
                        if (squaredDist < minSquaredDist) {
                            minSquaredDist = squaredDist;
                            nextIndexO = i;
                        }
                    }
                }
            }
            const size_t nextMolecule = nextIndexO - getStartIndexO();
            for (size_t i = 0; i < ring.getLength(); ++i) {
                if (nextMolecule == ring[i]) {
                    ringStart = i;
                    isRingGenerated = true;
                    goto done;
                }
            }
            ring.append(nextMolecule);
        done:;
        }

        Array<size_t> result(ring.getLength() - ringStart);
        for (size_t i = 0; i < result.getLength(); ++i)
            result[i] = ring[i + ringStart];
        return result;
    }
    /**
     * Reference:
     * [1] J. Chem. Phys. 118, 9291 (2003); https://doi.org/10.1063/1.1568337
     */
    template<Scalar T>
    auto IceGenerator<T>::makeRingMove(const Array<size_t>& ring, PositionMatrix& momentumMat) const -> CrystalCellType {
        const bool isInvalidRing = ring.getLength() == 0;
        if (isInvalidRing)
            return initialCell;

        PositionMatrix pos = initialCell.getPos();
        for (size_t bead = 0; bead < ring.getLength(); ++bead) {
            const size_t lastBead = bead == 0 ? ring.getLength() - 1 : bead - 1;
            const size_t nextBead = bead == ring.getLength() - 1 ? 0 : bead + 1;
            const size_t indexO = getStartIndexO() + ring[bead];
            const size_t indexLastO = getStartIndexO() + ring[lastBead];
            const size_t indexNextO = getStartIndexO() + ring[nextBead];
            const Vector3D<T> vecO2O1 = initialCell.minDistVector(indexO, indexLastO);
            const Vector3D<T> vecO2O3 = initialCell.minDistVector(indexO, indexNextO);

            const auto pair = findHydrogenInMolecule(ring[bead]);
            const size_t indexH1 = findHydrogenBetweenO(indexO, indexNextO);
            const size_t indexH2 = pair.first == indexH1 ? pair.second : pair.first;
            assert(indexH1 == pair.first || indexH1 == pair.second);

            const Vector3D<T> vecO2H2 = initialCell.minDistVector(indexO, indexH2);
            T angle;
            {
                const Vector3D<T> v1 = vecO2H2.cross(vecO2O1);
                const Vector3D<T> v2 = vecO2H2.cross(vecO2O3);
                angle = arccos(v1 * v2 / (v1.norm() * v2.norm()));
            }
            const Vector3D<T> vecO2H1 = initialCell.minDistVector(indexO, indexH1);
            const Vector3D<T> cross = vecO2O1.cross(vecO2O3);
            const T dot = cross * vecO2H2;
            if (dot.isNegative())
                angle = -angle;
            pos.row(indexH1) = pos.row(indexO) + rotate(angle, vecO2H2, vecO2H1);
            momentumMat.row(indexH1) = rotate(angle, vecO2H2, momentumMat.row(indexH1));
        }
        CrystalCellType result({initialCell.getLattice(), std::move(pos), CrystalCellType::Type::Cartesian}, initialCell.getAtomicNumbers());
        result.normalize();
        return result;
    }

    template<Scalar T>
    void IceGenerator<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        initialCell.swap(obj.initialCell);
        maxDistOO.swap(obj.maxDistOO);
        maxDistOH.swap(obj.maxDistOH);
        isHydrogenOccupied.swap(obj.isHydrogenOccupied);
        numHydrogenRequired.swap(obj.numHydrogenRequired);
    }

    template<Scalar T>
    void IceGenerator<T>::setInitialCell(CrystalCellType cell) {
        initialCell.swap(cell);
        if (initialCell.getType() == CrystalCellType::Type::Direct)
            initialCell.toCartesian();
        isHydrogenOccupied.resize(getEndIndexH());
        numHydrogenRequired.resize(getNumMolecule());
    }

    template<Scalar T>
    auto IceGenerator<T>::prepareRun() -> PositionMatrix {
        for (auto& elem : isHydrogenOccupied)
            elem = false;
        for (auto& elem : numHydrogenRequired)
            elem = 2U;
        return initialCell.getPos();
    }

    template<Scalar T>
    void IceGenerator<T>::searchDanglingH(PositionMatrix& pos) {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto pair = findHydrogenInMolecule(i);
            const auto hInRadius = findBondedH(i);
            const bool isDanglingH1 = std::find(hInRadius.begin(), hInRadius.end(), pair.first) == hInRadius.end();
            const bool isDanglingH2 = std::find(hInRadius.begin(), hInRadius.end(), pair.second) == hInRadius.end();
            if (isDanglingH1)
                fetchHydrogen(pos, i, pair.first);

            if (isDanglingH2)
                fetchHydrogen(pos, i, pair.second);
        }
    }

    template<Scalar T>
    Array<size_t> IceGenerator<T>::findOInRadius(size_t indexO, T radius) const {
        Array<size_t> result{};
        const T squaredRadiusO = square(radius);
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

    template<Scalar T>
    Array<size_t> IceGenerator<T>::findBondedH(size_t indexMolecule) const {
        const auto range = CrystalCellType::estimateRange(initialCell.getLattice(), maxDistOO);
        Array<size_t> result{};
        result.reserve(getNumMolecule() * 2);
        CrystalCellType::forCellInRange(range, initialCell.getLattice(), [this, indexMolecule, &result](Vector3D<T> delta) {
            const T squaredRadiusH = square(maxDistOO * 0.5);
            const T squaredRadiusO = square(maxDistOO);
            const size_t indexO = getStartIndexO() + indexMolecule;

            for (size_t i = getStartIndexO(); i < getEndIndexO(); ++i) {
                const Vector3D<T> otherO = initialCell.getPos().row(i) + delta;
                const T r2 = (initialCell.getPos().row(indexO) - otherO).squaredNorm();
                const bool isNotSelf = r2 > T(std::numeric_limits<T>::epsilon());
                const bool isInRange = r2 < squaredRadiusO;
                if (isNotSelf && isInRange) {
                    const Vector3D<T> middle = (otherO + initialCell.getPos().row(indexO)) * T(0.5);
                    for (size_t j = 0; j < getEndIndexH(); ++j) {
                        const bool isHInRange = initialCell.minDistVector(middle, j).squaredNorm() < squaredRadiusH;
                        const bool isOccupied = std::ranges::find(result, j) != result.end();
                        if (isHInRange && !isOccupied)
                            result.append(j);
                    }
                }
            }
        });
        return result;
    }

    template<Scalar T>
    Array<size_t> IceGenerator<T>::findFreeBondedHInRadius(size_t indexO) const {
        Array<size_t> result{};
        const auto hInRange = findBondedH(indexO);
        for (auto h : hInRange)
            if (!isHydrogenOccupied[h])
                result.append(h);
        return result;
    }

    template<Scalar T>
    std::pair<size_t, size_t> IceGenerator<T>::findHydrogenInMolecule(size_t indexMolecule) const {
        T minSquaredDist1 = std::numeric_limits<T>::max();
        T minSquaredDist2 = std::numeric_limits<T>::max();
        size_t indexH1 = 0;
        size_t indexH2 = 1;
        for (size_t i = 0; i < getEndIndexH(); ++i) {
            const T squaredDist = initialCell.minDistVector(i, getStartIndexO() + indexMolecule).squaredNorm();
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

    template<Scalar T>
    size_t IceGenerator<T>::findHydrogenBetweenO(size_t indexO1, size_t indexO2) const {
        T minSquaredDist = std::numeric_limits<T>::max();
        const Vector3D<T> delta = initialCell.minDistVector(indexO1, indexO2);
        const Vector3D<T> middle = initialCell.getPos().row(indexO1) + delta * T(0.5);
        size_t result = 0;
        for (size_t i = 0; i < getEndIndexH(); ++i) {
            const T squaredDist = initialCell.minDistVector(middle, i).squaredNorm();
            if (squaredDist < minSquaredDist) {
                minSquaredDist = squaredDist;
                result = i;
            }
        }
        return result;
    }

    template<Scalar T>
    template<RNG R>
    size_t IceGenerator<T>::makeRandEmptyO() const {
        std::uniform_int_distribution<size_t> dist(0, getNumMolecule() - 1);
        auto& gen = R::getInstance();
        size_t result = dist(gen);
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto hInRange = findBondedH(result);
            const size_t numFreeH = countFreeH(hInRange);
            const bool isOxygenEmpty = numHydrogenRequired[result] != 0;
            if (isOxygenEmpty && numFreeH != 0)
                break;
            result = (result + 1U) % getNumMolecule();
        }
        assert(numHydrogenRequired[result] != 0);
        return result;
    }

    template<Scalar T>
    template<RNG R>
    size_t IceGenerator<T>::makeRandFreeH(size_t indexO) const {
        assert(indexO < getNumMolecule());
        const auto hInRange = findBondedH(indexO);
        const size_t randLogicIndex = [&]() {
            const size_t numFreeH = countFreeH(hInRange);
            assert(numFreeH > numHydrogenRequired[indexO]);
            std::uniform_int_distribution<size_t> dist(0, numFreeH - 1);
            return dist(R::getInstance());
        }();

        size_t logicIndex = 0;
        for (size_t physicalIndex = 0; physicalIndex < hInRange.getLength(); ++physicalIndex) {
            const bool isFree = isHydrogenOccupied[hInRange[physicalIndex]] == false;
            if (isFree && logicIndex == randLogicIndex)
                return hInRange[physicalIndex];
            logicIndex += isFree;
        }
        unreachable();
    }

    template<Scalar T>
    void IceGenerator<T>::fetchHydrogen(PositionMatrix& pos, size_t indexO, size_t indexH) {
        assert(!isHydrogenOccupied[indexH]);
        assert(numHydrogenRequired[indexO] > 0);
        
        const size_t shiftIndexO = getStartIndexO() + indexO;
        Vector3D<T> delta = initialCell.minDistVector(shiftIndexO, indexH);
        delta.toUnit();
        delta *= T(BondLengthOH);
        auto row = pos.row(indexO * 2U + (2U - numHydrogenRequired[indexO]));
        row = initialCell.getPos().row(shiftIndexO) + delta;

        isHydrogenOccupied[indexH] = true;
        numHydrogenRequired[indexO] -= 1;
    }

    template<Scalar T>
    size_t IceGenerator<T>::countFreeH(const Array<size_t>& hIndexes) const {
        size_t numFreeH = 0;
        for (auto h : hIndexes)
            numFreeH += !isHydrogenOccupied[h];
        return numFreeH;
    }

    template<Scalar T>
    void IceGenerator<T>::searchForPairs(PositionMatrix& pos) {
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

    template<Scalar T>
    size_t IceGenerator<T>::getIndexToPair() const {
        for (size_t indexO = 0; indexO < getNumMolecule(); ++indexO) {
            const auto hInRange = findBondedH(indexO);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH > 0 && numFreeH == numHydrogenRequired[indexO])
                return indexO;
        }
        return getNumMolecule();
    }

    template<Scalar T>
    template<RNG R>
    void IceGenerator<T>::randUninitializedH(PositionMatrix& pos) {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            while (numHydrogenRequired[i] != 0) {
                auto row = pos.row(i * 2U + (2U - numHydrogenRequired[i]));
                row = initialCell.getPos().row(i + getStartIndexO()) + randUnitVector<R>() * T(BondLengthOH);
                numHydrogenRequired[i] -= 1;
            }
        }
    }

    template<Scalar T>
    void IceGenerator<T>::exhaustImpl(size_t stackDepth, const PositionMatrix& pos, Array<CrystalCellType>& result) {
        const bool recursionStop = stackDepth == getNumMolecule();
        if (recursionStop) {
            CrystalCellType cell({initialCell.getLattice(), pos, CrystalCellType::Type::Cartesian}, initialCell.getAtomicNumbers());
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

    template<Scalar T>
    bool IceGenerator<T>::isFinished() const noexcept {
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            const auto hInRange = findBondedH(i);
            const size_t numFreeH = countFreeH(hInRange);
            if (numFreeH != 0 && numHydrogenRequired[i] != 0)
                return false;
        }
        return true;
    }

    template<Scalar T>
    Vector3D<T> IceGenerator<T>::rotate(T angle, Vector3D<T> axis, Vector3D<T> target) {
        const T factor1 = cos(angle);

        const Vector3D<T> parallel = (target * axis / axis.squaredNorm()) * axis;
        const Vector3D<T> verticle = target - parallel;
        const Vector3D<T> verticle2 = parallel.cross(verticle);
        T factor2 = sin(angle) / parallel.norm();
        return Vector3D<T>(parallel + factor1 * verticle + factor2 * verticle2);
    }

    template<Scalar T>
    template<RNG R>
    Vector3D<T> IceGenerator<T>::randUnitVector() {
        const T theta(T::template random_uniform<R>() * M_PI);
        const T phi(T::template random_uniform<R>() * M_PI * 2);
        Vector3D<T> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
        return result;
    }
}
