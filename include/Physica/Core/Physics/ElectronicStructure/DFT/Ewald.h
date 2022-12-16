/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"
#include "Physica/Utils/TestHelper.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     * [3] Toukmaji A Y, Board J A. Ewald summation techniques in perspective: a survey[J]. Computer Physics Communications, 1996, 95(2-3):73-92.
     */
    template<class ScalarType, class PosScalarType = ScalarType>
    class Ewald {
        constexpr static unsigned int Dim = 3;
        using LatticeMatrix = typename PeriodicCell<PosScalarType, Dim>::LatticeMatrix;
        using PositionMatrix = typename PeriodicCell<PosScalarType, Dim>::PositionMatrix;
        using SearchRangeType = typename PeriodicCell<PosScalarType, Dim>::SearchRangeType;
        using CellListType = CellList<ScalarType, PosScalarType>;
        using Index3D = typename CellListType::Index3D;
        using Vector3D = Vector<PosScalarType, Dim>;
        constexpr static size_t ErfcTableSize = 4096 + 512;
        constexpr static double ErfcTableStep = 0.001;
        constexpr static double SumPrec = ErfcTableSize * ErfcTableStep; // Refer to [2]
    private:
        LatticeMatrix lattice;
        ReciprocalCell<PosScalarType> repCell;
        Vector<ScalarType> charges;
        ScalarType inv_volume;
        ScalarType integralLimit;
        ScalarType selfEnergy;
        PosScalarType rSpaceCutoff;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
        Vector<ScalarType> erfc_table;
        ScalarType erfcStep;
        ScalarType repErfcStep;
        ScalarType doubleSquareStep;
        ScalarType maxErfcX;
        ScalarType halfErfcStep;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_);
        Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<ScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const PositionMatrix& pos) const;
        [[nodiscard]] ScalarType potentialEnergy(const PositionMatrix& pos) const;
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
    private:
        void initIntegralLimit(PosScalarType volume);
        void makeTables();
        inline ScalarType calcFromTable(ScalarType x) const;
        inline ScalarType calcFromTable_diff(ScalarType x) const;
        Vector<ScalarType> rSpaceForce(const PositionMatrix& pos) const;

        friend class ::Physica::Test;
    };

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_)
            : Ewald(lattice_, ReciprocalCell(lattice_), std::move(charges_)) {}

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<ScalarType> charges_)
            : lattice(std::move(lattice_))
            , repCell(std::move(repCell_))
            , charges(std::move(charges_))
            , erfc_table(ErfcTableSize) {
        const PosScalarType volume = PeriodicCell<PosScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(ScalarType(volume));
        initIntegralLimit(volume);
        rSpaceCutoff = ScalarType(SumPrec) / integralLimit;
        rSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(lattice, rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(repCell.getLattice(), PosScalarType(ScalarType(SumPrec * 2) * integralLimit));
        selfEnergy = square(charges).sum() * integralLimit / sqrt(ScalarType(M_PI))
                   + square(charges.sum()) * ScalarType(M_PI) / (ScalarType::Two() * square(integralLimit)) * inv_volume;
        makeTables();
    }

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>& Ewald<ScalarType, PosScalarType>::operator=(Ewald ewald) noexcept {
        swap(ewald);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType, PosScalarType>::force(const PositionMatrix& pos) const {
        Vector<ScalarType> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            const size_t numParticle = getNumParticle();
            Vector<ScalarType> kSpaceSum(numParticle * Dim, 0);
            const ScalarType factor1 = reciprocal(square(ScalarType::Two() * integralLimit));
            Vector<ScalarType> dots(numParticle);
            Vector<ScalarType> sin_vec(numParticle);
            Vector<ScalarType> cos_vec(numParticle);
            PeriodicCell<PosScalarType, Dim>::forReducedCellInRange(kSpaceSumRange, repCell.getLattice(),
                [this, numParticle, factor1, &dots, &sin_vec, &cos_vec, &pos, &kSpaceSum](Vector3D delta) {
                    const ScalarType squaredNorm = ScalarType(delta.squaredNorm());
                    const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                    if (isNotGammaPoint) {
                        dots = pos * delta;
                        sincos(dots, sin_vec, cos_vec);
                        const ScalarType sum_cos = cos_vec * charges;
                        const ScalarType sum_sin = sin_vec * charges;
                        const ScalarType factor2 = reciprocal(squaredNorm * exp(squaredNorm * factor1));
                        for (size_t i = 0; i < numParticle; ++i) {
                            auto force_i = kSpaceSum.segment(i * Dim, (i + 1) * Dim);
                            const ScalarType charge = charges[i];
                            force_i += ((sin_vec[i] * sum_cos - cos_vec[i] * sum_sin) * (factor2 * charge)) * delta;
                        }
                    }
                });
            kSpaceSum *= ScalarType(8 * M_PI) * inv_volume;
            result = std::move(kSpaceSum);
        });
        Vector<ScalarType> rSpaceSum = rSpaceForce(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class ScalarType, class PosScalarType>
    ScalarType Ewald<ScalarType, PosScalarType>::potentialEnergy(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        ScalarType rSpaceSum = 0;
        PeriodicCell<PosScalarType, Dim>::forCellInRange(rSpaceSumRange, lattice, [this, numParticle, &pos, &rSpaceSum](Vector3D delta) {
                ScalarType sum = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    const Vector3D pos_i = pos.row(i).asVector() + delta;
                    const ScalarType charge_i = charges[i];
                    for (size_t j = i; j < numParticle; ++j) {
                        const Vector3D pos_ij = pos_i - pos.row(j).asVector();
                        const ScalarType r2 = pos_ij.squaredNorm();
                        const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                        if (isNotSelf) {
                            const ScalarType r = sqrt(r2);
                            const ScalarType temp = calcFromTable(r) * (charge_i * charges[j]);
                            sum += temp * ScalarType(i == j ? 1 : 2);
                        }
                    }
                }
                rSpaceSum += sum;
            });
        
        ScalarType kSpaceSum = 0;
        const ScalarType factor = reciprocal(square(ScalarType::Two() * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forCellInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor, &pos, &kSpaceSum](Vector3D delta) {
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                if (isNotGammaPoint) {
                    ScalarType sum_cos = 0;
                    ScalarType sum_sin = 0;
                    for (size_t i = 0; i < numParticle; ++i) {
                        const ScalarType charge = charges[i];
                        const ScalarType dot = delta * pos.row(i).asVector();
                        ScalarType cos_temp, sin_temp;
                        sincos(dot, sin_temp, cos_temp);
                        sum_cos += charge * cos_temp;
                        sum_sin += charge * sin_temp;
                    }
                    kSpaceSum += (square(sum_cos) + square(sum_sin)) / (squaredNorm * exp(squaredNorm * factor));
                }
            });
        kSpaceSum *= ScalarType(4 * M_PI) * inv_volume;
        return (rSpaceSum + kSpaceSum) * 0.5 - selfEnergy;
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::swap(Ewald& ewald) noexcept {
        lattice.swap(ewald.lattice);
        repCell.swap(ewald.repCell);
        charges.swap(ewald.charges);
        inv_volume.swap(ewald.inv_volume);
        integralLimit.swap(ewald.integralLimit);
        selfEnergy.swap(ewald.selfEnergy);
        rSpaceCutoff.swap(ewald.rSpaceCutoff);
        rSpaceSumRange.swap(ewald.rSpaceSumRange);
        kSpaceSumRange.swap(ewald.kSpaceSumRange);
        erfc_table.swap(ewald.erfc_table);
        erfcStep.swap(ewald.erfcStep);
        repErfcStep.swap(ewald.repErfcStep);
        doubleSquareStep.swap(ewald.doubleSquareStep);
        maxErfcX.swap(ewald.maxErfcX);
        halfErfcStep.swap(ewald.halfErfcStep);
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::initIntegralLimit(PosScalarType volume) {
        const ScalarType averageCellSize = cbrt(ScalarType(volume));
        const ScalarType estimate = sqrt(cbrt(ScalarType(getNumParticle())) * M_PI) / averageCellSize;

        constexpr double factor = 2 * M_PI;
        const auto& repLatt = repCell.getLattice();
        const ScalarType heightX = reciprocal(repLatt.row(0).norm()) * factor;
        const ScalarType heightY = reciprocal(repLatt.row(1).norm()) * factor;
        const ScalarType heightZ = reciprocal(repLatt.row(2).norm()) * factor;
        constexpr double factor1 = 0.99; //A rule of thumb
        const ScalarType maxRSpaceCutoff = std::min(heightX, std::min(heightY, heightZ)) * factor1;
        const ScalarType minLimit = ScalarType(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(estimate, minLimit);
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::makeTables() {
        for (size_t i = 0; i < erfc_table.getLength(); ++i) {
            ScalarType x = ScalarType(i * ErfcTableStep);
            erfc_table[i] = erfc(x) / x * integralLimit;
        }
        erfcStep = ScalarType(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        doubleSquareStep = square(erfcStep) * 2;
        maxErfcX = (ScalarType(ErfcTableSize - 1) - 0.5) * erfcStep;
        halfErfcStep = erfcStep * 0.5;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType Ewald<ScalarType, PosScalarType>::calcFromTable(ScalarType x) const {
        if (x > maxErfcX)
            return 0;
        const bool safeToRound = x > halfErfcStep;
        if (safeToRound) {
            const ScalarType temp = x * repErfcStep + 0.5;
            const size_t index = double(temp);
            const ScalarType x2 = erfcStep * floor(temp);
            auto y = erfc_table.segment(index - 1, index + 2);
            return Internal::quadraticInterpolate(x2 - erfcStep, x2, x2 + erfcStep, y[0], y[1], y[2], x); //Optimize: make use of x1, x2, x3 are equal distance
        }
        return 1;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType Ewald<ScalarType, PosScalarType>::calcFromTable_diff(ScalarType x) const {
        if (x > maxErfcX)
            return 0;
        const bool safeToRound = x > halfErfcStep;
        if (safeToRound) {
            const ScalarType temp = x * repErfcStep + 0.5;
            const size_t index = double(temp);
            const ScalarType x2 = erfcStep * floor(temp);
            auto y = erfc_table.segment(index - 1, index + 2);
            const ScalarType factor = doubleSquareStep * x;
            return Internal::quadraticInterpolate_diff1(factor, erfcStep, x2, y[0], y[1], y[2], x);
        }
        return 1;
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> Ewald<ScalarType, PosScalarType>::rSpaceForce(const PositionMatrix& pos) const {
        const auto numParticle = getNumParticle();
        Vector<ScalarType> rSpaceSum(numParticle * Dim, 0);
        const CellListType cellList(lattice, pos, rSpaceCutoff);
        for (size_t i = 0; i < numParticle; ++i) {
            const Index3D center = cellList.getAtomCellMap()[i];
            cellList.forNeighInRange(center, [this, i, pos, &rSpaceSum, &cellList](Vector3D translate, Index3D neigh) {
                const Vector3D from = pos.row(i) - translate;
                Vector3D delta;
                Vector<ScalarType, Dim> sum(Dim, 0);
                for (size_t j : cellList(neigh)) {
                    const auto to = pos.row(j);
                    delta = from - to;
                    const ScalarType charge = charges[j];
                    const ScalarType temp = charge * calcFromTable_diff(delta.norm());
                    sum += temp * delta;
                }
                auto f = rSpaceSum.segment(i * Dim, (i + 1) * Dim);
                f += charges[i] * sum;
            });
        }
        return rSpaceSum;
    }
}
