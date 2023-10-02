/*
 * Copyright 2021-2023 WeiBo He.
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

#include <iostream>
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Utils/TestHelper.h"
#include "PairModel.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType> class Ewald;

    namespace Internal {
        template<class T, class U>
        class Traits<Ewald<T, U>> {
        public:
            using ScalarType = T;
            using PosScalarType = U;
            constexpr static bool IsPotDependOnAtomIndex = true;
        };
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     * [3] Toukmaji A Y, Board J A. Ewald summation techniques in perspective: a survey[J]. Computer Physics Communications, 1996, 95(2-3):73-92.
     */
    template<class ScalarType, class PosScalarType>
    class Ewald : public PairModel<Ewald<ScalarType, PosScalarType>> {
        constexpr static unsigned int Dim = 3;
        using Base = PairModel<Ewald<ScalarType, PosScalarType>>;
        using PlainScalar = typename ScalarType::PlainScalar;
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
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
        Vector<PlainScalar> erfc_table;
        ScalarType erfcStep;
        ScalarType repErfcStep;
        ScalarType repDoubleSquareStep;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const PositionMatrix& pos) const;
        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] inline Vector<ScalarType> force_short(const PositionMatrix& pos) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const PositionMatrix& pos) const;

        [[nodiscard]] ScalarType potentialEnergy(const PositionMatrix& pos) const;
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] PosScalarType getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] PosScalarType getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        /* Setters */
        void setLattice(LatticeMatrix lattice_);
    protected:
        inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    private:
        void initIntegralLimit(PosScalarType volume);
        void makeTables();

        using Base::getCutoff;
        using Base::getSquaredCutoff;

        friend class ::Physica::Test;
        friend class PairModel<Ewald<ScalarType, PosScalarType>>;
    };

    template<class ScalarType, class PosScalarType>
    Ewald<ScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_)
            : charges(std::move(charges_))
            , erfc_table(ErfcTableSize + 1) {
        setLattice(std::move(lattice_));
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
            result = force_long<SequentialExecutor>(pos);
        });
        const Vector<ScalarType> rSpaceSum = force_short<SequentialExecutor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool IsSmallCell> 
    inline Vector<ScalarType> Ewald<ScalarType, PosScalarType>::force_short(const PositionMatrix& pos) const {
        static_assert(std::is_same<Executor, SequentialExecutor>::value, "[Error]: Parallelization is not implemented");
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        const Vector<ScalarType> rSpaceSum = Base::template force<SequentialExecutor, false>(lattice, pos);
        return rSpaceSum;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType, PosScalarType>::force_long(const PositionMatrix& pos) const {
        static_assert(std::is_same<Executor, SequentialExecutor>::value, "[Error]: Parallelization is not implemented");
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> kSpaceSum(numParticle * Dim, 0);
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        const ScalarType factor1 = reciprocal(square(PlainScalar(2) * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forReducedCellInRange(kSpaceSumRange, repCell.getLattice(), // Reduce cell using time reversal symmetry
            [this, numParticle, factor1, &dots, &sin_vec, &cos_vec, &pos, &kSpaceSum](Vector3D delta) {
                const ScalarType squaredNorm = ScalarType(delta.squaredNorm());
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const ScalarType sum_cos = cos_vec * charges;
                    const ScalarType sum_sin = sin_vec * charges;
                    const ScalarType factor2 = reciprocal(squaredNorm * exp(squaredNorm * factor1));
                    for (size_t i = 0; i < numParticle; ++i) {
                        auto force_i = kSpaceSum.template segment<3>(i * Dim, (i + 1) * Dim);
                        const ScalarType temp = (sin_vec[i] * sum_cos - cos_vec[i] * sum_sin) * (factor2 * charges[i]);
                        force_i[0] += temp * delta[0];
                        force_i[1] += temp * delta[1];
                        force_i[2] += temp * delta[2];
                    }
                }
            });
        kSpaceSum *= ScalarType(8 * M_PI) * inv_volume; // 8 Pi because time reversal symmetry
        return kSpaceSum;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class ScalarType, class PosScalarType>
    ScalarType Ewald<ScalarType, PosScalarType>::potentialEnergy(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        const ScalarType rSpaceSum = Base::potentialEnergy(lattice, pos);
        ScalarType kSpaceSum = 0;
        const ScalarType factor = reciprocal(square(PlainScalar(2) * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forCellInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor, &pos, &kSpaceSum](Vector3D delta) {
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
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
        kSpaceSum *= PlainScalar(4 * M_PI) * inv_volume;
        return (rSpaceSum + kSpaceSum) * PlainScalar(0.5) - selfEnergy;
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::swap(Ewald& ewald) noexcept {
        assert(this != &ewald && "[Error]: Self swap is likely a bug");
        Base::swap(ewald);
        lattice.swap(ewald.lattice);
        repCell.swap(ewald.repCell);
        charges.swap(ewald.charges);
        inv_volume.swap(ewald.inv_volume);
        integralLimit.swap(ewald.integralLimit);
        selfEnergy.swap(ewald.selfEnergy);
        rSpaceSumRange.swap(ewald.rSpaceSumRange);
        kSpaceSumRange.swap(ewald.kSpaceSumRange);
        erfc_table.swap(ewald.erfc_table);
        erfcStep.swap(ewald.erfcStep);
        repErfcStep.swap(ewald.repErfcStep);
        repDoubleSquareStep.swap(ewald.repDoubleSquareStep);
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::setLattice(LatticeMatrix lattice_) {
        lattice = std::move(lattice_);
        repCell = ReciprocalCell(lattice);
        const PosScalarType volume = PeriodicCell<PosScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(ScalarType(volume));
        initIntegralLimit(volume);
        const PlainScalar rSpaceCutoff = PlainScalar(SumPrec) / integralLimit.getValue();
        rSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(lattice, rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(repCell.getLattice(), PlainScalar(SumPrec * 2) * integralLimit.getValue());
        selfEnergy = square(charges).sum() * integralLimit / sqrt(PlainScalar(M_PI))
                   + square(charges.sum()) * PlainScalar(M_PI) / (PlainScalar(2) * square(integralLimit)) * inv_volume;
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType Ewald<ScalarType, PosScalarType>::force_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = double(temp);
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        return charges[i] * charges[j] * Internal::quadraticInterpolate_diff1<ScalarType>(repDoubleSquareStep, erfcStep, x1, y[0], y[1], y[2], r);
    }
    /**
     * Optimize: make use of x1, x2, x3 are equal distance
     */
    template<class ScalarType, class PosScalarType>
    inline ScalarType Ewald<ScalarType, PosScalarType>::pot_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = double(temp);
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const ScalarType interp = Internal::quadraticInterpolate<ScalarType>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return ScalarType(i == j ? 1 : 2) * charges[i] * charges[j] * interp;
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::initIntegralLimit(PosScalarType volume) {
        const ScalarType averageCellSize = cbrt(ScalarType(volume));
        const ScalarType estimate = sqrt(cbrt(ScalarType(getNumParticle())) * ScalarType(M_PI)) / averageCellSize;

        const auto& repLatt = repCell.getLattice();
        const ScalarType heightX_2Pi = reciprocal(repLatt.row(0).norm());
        const ScalarType heightY_2Pi = reciprocal(repLatt.row(1).norm());
        const ScalarType heightZ_2Pi = reciprocal(repLatt.row(2).norm());
        constexpr double factor1 = 2 * M_PI * (1 - std::numeric_limits<ScalarType>::epsilon()); //To avoid rSpaceCutoff larger than max value
        const ScalarType maxRSpaceCutoff = std::min(heightX_2Pi, std::min(heightY_2Pi, heightZ_2Pi)) * PlainScalar(factor1);
        const ScalarType minLimit = ScalarType(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(estimate, minLimit);
    }

    template<class ScalarType, class PosScalarType>
    void Ewald<ScalarType, PosScalarType>::makeTables() {
        for (size_t i = 1; i < erfc_table.getLength(); ++i) {
            const auto x = PlainScalar((i - 1) * ErfcTableStep);
            erfc_table[i] = erfc(x) / x * integralLimit.getValue();
        }
        erfc_table[0] = erfc_table[1] = erfc_table[2]; // Smooth out divergent erfc(0) / 0 
        erfcStep = PlainScalar(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * PlainScalar(2));
    }
}
