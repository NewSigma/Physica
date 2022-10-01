/*
 * Copyright 2021 WeiBo He.
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
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Grid3D.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     * [3] VASP (www.vasp.at)
     * [4] Toukmaji A Y, Board J A. Ewald summation techniques in perspective: a survey[J]. Computer Physics Communications, 1996, 95(2-3):73-92.
     */
    template<class ScalarType, class PosScalarType = ScalarType>
    class Ewald {
        constexpr static unsigned int Dim = 3;
        using LatticeMatrix = typename PeriodicCell<PosScalarType, Dim>::LatticeMatrix;
        using PositionMatrix = typename PeriodicCell<PosScalarType, Dim>::PositionMatrix;
        using SearchRangeType = typename PeriodicCell<PosScalarType, Dim>::SearchRangeType;
        using Vector3D = Vector<PosScalarType, Dim>;
        constexpr static size_t ErfcTableSize = 4096 + 512;
        constexpr static double ErfcTableStep = 0.001;
        constexpr static size_t ExpTableSize = 8192 + 64;
        constexpr static double ExpTableStep = 0.001;
        constexpr static double SumPrec = 4; // Refer to [2]
    private:
        LatticeMatrix lattice;
        ReciprocalCell<PosScalarType> repCell;
        Vector<ScalarType> charges;
        ScalarType inv_volume;
        ScalarType integralLimit;
        ScalarType selfEnergy;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
        Vector<ScalarType> erfc_table;
        Vector<ScalarType> exp_table;
        ScalarType erfcStep;
        ScalarType repErfcStep;
        ScalarType expStep;
        ScalarType repExpStep;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_);
        Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<ScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        [[nodiscard]] Vector<ScalarType> force(const PositionMatrix& pos) const;
        [[nodiscard]] ScalarType potentialEnergy(const PositionMatrix& pos) const;
        /* Operations */
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
    private:
        void makeTables();
        static inline ScalarType calcFromTable(
                const Vector<ScalarType>& table,
                ScalarType step,
                ScalarType repStep,
                ScalarType x);
    };

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_)
            : Ewald(lattice_, ReciprocalCell(lattice_), std::move(charges_)) {}

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<ScalarType> charges_)
            : lattice(std::move(lattice_))
            , repCell(std::move(repCell_))
            , charges(std::move(charges_))
            , erfc_table(ErfcTableSize)
            , exp_table(ExpTableSize) {
        const PosScalarType volume = PeriodicCell<PosScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(ScalarType(volume));
        {
            const ScalarType averageCellSize = cbrt(ScalarType(volume));
            integralLimit = sqrt(ScalarType(M_PI)) / averageCellSize;
            rSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(lattice, PosScalarType(ScalarType(SumPrec) / integralLimit));
            kSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(repCell.getLattice(), PosScalarType(ScalarType(SumPrec * 2) * integralLimit));
        }
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
    Vector<ScalarType> Ewald<ScalarType, PosScalarType>::force(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> result(numParticle * Dim, 0);
        
        const ScalarType factor1 = reciprocal(square(ScalarType::Two() * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forCellInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor1, &pos, &result](Vector3D delta) {
                const ScalarType squaredNorm = ScalarType(delta.squaredNorm());
                const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                if (isNotGammaPoint) {
                    const ScalarType factor2 = reciprocal(squaredNorm * exp(squaredNorm * factor1));
                    ScalarType sum_cos = 0;
                    ScalarType sum_sin = 0;
                    for (size_t i = 0; i < numParticle; ++i) {
                        const ScalarType charge = charges[i];
                        const ScalarType dot = ScalarType(delta * pos.row(i).asVector());
                        ScalarType cos_temp, sin_temp;
                        sincos(dot, sin_temp, cos_temp);
                        sum_cos += charge * cos_temp;
                        sum_sin += charge * sin_temp;
                    }
                    for (size_t i = 0; i < numParticle; ++i) {
                        auto force_i = result.segment(i * Dim, (i + 1) * Dim);
                        const ScalarType charge = charges[i];
                        const ScalarType dot = ScalarType(delta * pos.row(i).asVector());
                        ScalarType cos_temp, sin_temp;
                        sincos(dot, sin_temp, cos_temp);
                        force_i += ((sin_temp * sum_cos - cos_temp * sum_sin) * (factor2 * charge)) * delta;
                    }
                }
            });
        result *= ScalarType(4 * M_PI) * inv_volume;

        PeriodicCell<PosScalarType, Dim>::forCellInRange(rSpaceSumRange, lattice, [this, numParticle, &pos, &result](Vector3D delta) {
                for (size_t i = 0; i < numParticle; ++i) {
                    const Vector3D pos_i = pos.row(i).asVector() + delta;
                    Vector<ScalarType, Dim> sum(Dim, 0);
                    for (size_t j = 0; j < numParticle; ++j) {
                        if (j == i)
                            continue;
                        const Vector3D pos_ij = pos_i - pos.row(j).asVector();
                        const ScalarType r2 = ScalarType(pos_ij.squaredNorm());
                        const ScalarType r = sqrt(r2);
                        const ScalarType temp1 = calcFromTable(exp_table, expStep, repExpStep, r2);
                        const ScalarType temp2 = calcFromTable(erfc_table, erfcStep, repErfcStep, r);
                        const ScalarType temp3 = (temp1 + temp2) * (charges[j] / r2);
                        sum += temp3 * pos_ij;
                    }
                    auto force_i = result.segment(i * Dim, (i + 1) * Dim);
                    force_i += charges[i] * sum;
                }
            });
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
                            const ScalarType temp = calcFromTable(erfc_table, erfcStep, repErfcStep, r) * (charge_i * charges[j]);
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
        rSpaceSumRange.swap(ewald.rSpaceSumRange);
        kSpaceSumRange.swap(ewald.kSpaceSumRange);
        erfc_table.swap(ewald.erfc_table);
        exp_table.swap(ewald.exp_table);
        erfcStep.swap(ewald.erfcStep);
        repErfcStep.swap(ewald.repErfcStep);
        expStep.swap(ewald.expStep);
        repExpStep.swap(ewald.repExpStep);
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::makeTables() {
        for (size_t i = 0; i < erfc_table.getLength(); ++i) {
            ScalarType x = ScalarType(i * ErfcTableStep);
            erfc_table[i] = erfc(x) / x * integralLimit;
        }
        erfcStep = ScalarType(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);

        for (size_t i = 0; i < exp_table.getLength(); ++i) {
            ScalarType x = ScalarType(i * ExpTableStep);
            exp_table[i] = exp(-x) * (ScalarType(2) / sqrt(ScalarType(M_PI)) * integralLimit);
        }
        expStep = ScalarType(ExpTableStep) / square(integralLimit);
        repExpStep = reciprocal(expStep);
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType Ewald<ScalarType, PosScalarType>::calcFromTable(
            const Vector<ScalarType>& table,
            ScalarType step,
            ScalarType repStep,
            ScalarType x) {
        const size_t index = double(x * repStep + 0.5);
        if (index == 0)
            return 1;
        if (index + 1 >= table.getLength())
            return 0;
        ScalarType x2 = step * index;
        auto y = table.segment(index - 1, index + 2);
        return Internal::quadraticInterpolate(x2 - step, x2, x2 + step, y[0], y[1], y[2], x);
    }
}
