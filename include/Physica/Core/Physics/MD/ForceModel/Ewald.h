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
    template<class ScalarType> class Ewald;

    namespace Internal {
        template<class T>
        class Traits<Ewald<T>> {
        public:
            using ScalarType = T;
            constexpr static bool IsPotDependOnAtomIndex = true;
        };
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     * [3] Toukmaji A Y, Board J A. Ewald summation techniques in perspective: a survey[J]. Computer Physics Communications, 1996, 95(2-3):73-92.
     */
    template<class ScalarType>
    class Ewald : public PairModel<Ewald<ScalarType>> {
        constexpr static unsigned int Dim = 3;
        using Base = PairModel<Ewald<ScalarType>>;
        using ComplexType = ComplexScalar<ScalarType>;
        using PlainScalar = typename ScalarType::PlainScalar;
        using LatticeMatrix = typename PeriodicCell<ScalarType, Dim>::LatticeMatrix;
        using PositionMatrix = typename PeriodicCell<ScalarType, Dim>::PositionMatrix;
        using SearchRangeType = typename PeriodicCell<ScalarType, Dim>::SearchRangeType;
        using CellListType = CellList<ScalarType>;
        using Index3D = typename CellListType::Index3D;
        using Vector3D = Vector<ScalarType, Dim>;
        constexpr static size_t ErfcTableSize = 4096 + 512 + 2;
        constexpr static double ErfcTableStep = 0.001;
        constexpr static double SumPrec = (ErfcTableSize - 2) * ErfcTableStep; // Referenced from [2], minus 2 to avoid overflow
    private:
        LatticeMatrix lattice;
        ReciprocalCell<ScalarType> repCell;
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
        [[nodiscard]] ComplexType forceConst(const PositionMatrix& pos, Vector3D qPoint, size_t dof1, size_t dof2) const;
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] ScalarType getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] ScalarType getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        /* Setters */
        void setLattice(LatticeMatrix lattice_);
    protected:
        inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    private:
        /* Operations */
        void initIntegralLimit(ScalarType volume);
        void makeTables();
        [[nodiscard]] ComplexType kSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const;
        [[nodiscard]] ComplexType rSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const;
        [[nodiscard]] ScalarType rSpaceForceConstImpl1(ScalarType r) const;
        [[nodiscard]] ScalarType rSpaceForceConstImpl2(ScalarType r) const;
        /* Getters */
        using Base::getCutoff;
        using Base::getSquaredCutoff;
        /* Friends */
        friend class ::Physica::Test;
        friend class PairModel<Ewald<ScalarType>>;
    };

    template<class ScalarType>
    Ewald<ScalarType>::Ewald(LatticeMatrix lattice_, Vector<ScalarType> charges_)
            : charges(std::move(charges_))
            , erfc_table(ErfcTableSize + 1) {
        setLattice(std::move(lattice_));
    }

    template<class ScalarType>
    Ewald<ScalarType>& Ewald<ScalarType>::operator=(Ewald ewald) noexcept {
        swap(ewald);
        return *this;
    }

    template<class ScalarType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType>::force(const PositionMatrix& pos) const {
        Vector<ScalarType> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            result = force_long<SequentialExecutor>(pos);
        });
        const Vector<ScalarType> rSpaceSum = force_short<SequentialExecutor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }

    template<class ScalarType>
    template<class Executor, bool IsSmallCell> 
    inline Vector<ScalarType> Ewald<ScalarType>::force_short(const PositionMatrix& pos) const {
        static_assert(std::is_same<Executor, SequentialExecutor>::value, "[Error]: Parallelization is not implemented");
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        const Vector<ScalarType> rSpaceSum = Base::template force<SequentialExecutor, false>(lattice, pos);
        return rSpaceSum;
    }

    template<class ScalarType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType>::force_long(const PositionMatrix& pos) const {
        static_assert(std::is_same<Executor, SequentialExecutor>::value, "[Error]: Parallelization is not implemented");
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> kSpaceSum(numParticle * Dim, 0);
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        const ScalarType factor1 = reciprocal(square(PlainScalar(2) * integralLimit));
        PeriodicCell<ScalarType, Dim>::forReducedCellInRange(kSpaceSumRange, repCell.getLattice(), // Reduce cell using time reversal symmetry
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
        if constexpr (ScalarType::isReverseDiff)
            kSpaceSum.makeContinuous();
        kSpaceSum *= ScalarType(8 * M_PI) * inv_volume; // 8 Pi because time reversal symmetry
        return kSpaceSum;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class ScalarType>
    ScalarType Ewald<ScalarType>::potentialEnergy(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        const ScalarType rSpaceSum = Base::potentialEnergy(lattice, pos);
        ScalarType kSpaceSum = 0;
        const ScalarType factor = reciprocal(square(PlainScalar(2) * integralLimit));
        PeriodicCell<ScalarType, Dim>::forCellInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor, &pos, &kSpaceSum](Vector3D delta) {
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
        return kSpaceSum * PlainScalar(0.5) + rSpaceSum - selfEnergy;
    }
    /**
     * Reference:
     * [1] Rev. Mod. Phys. 73, 515; https://doi.org/10.1103/RevModPhys.73.515
     */
    template<class ScalarType>
    typename Ewald<ScalarType>::ComplexType
    Ewald<ScalarType>::forceConst(const PositionMatrix& pos, Vector3D qPoint, size_t dof1, size_t dof2) const {
        const Vector3D waveQ = repCell.getLattice().transpose() * qPoint;
        const ComplexType kSpaceSum = kSpaceForceConst(pos, waveQ, dof1, dof2);
        const ComplexType rSpaceSum = rSpaceForceConst(pos, waveQ, dof1, dof2);
        return kSpaceSum + rSpaceSum;
    }

    template<class ScalarType>
    void Ewald<ScalarType>::swap(Ewald& ewald) noexcept {
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

    template<class ScalarType>
    void Ewald<ScalarType>::setLattice(LatticeMatrix lattice_) {
        lattice = std::move(lattice_);
        repCell = ReciprocalCell(lattice);
        const ScalarType volume = PeriodicCell<ScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(ScalarType(volume));
        initIntegralLimit(volume);
        const PlainScalar rSpaceCutoff = PlainScalar(SumPrec) / integralLimit.getValue();
        rSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(lattice, rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(repCell.getLattice(), PlainScalar(SumPrec * 2) * integralLimit.getValue());
        selfEnergy = square(charges).sum() * integralLimit / sqrt(PlainScalar(M_PI))
                   + square(charges.sum()) * PlainScalar(M_PI) / (PlainScalar(2) * square(integralLimit)) * inv_volume;
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<class ScalarType>
    inline ScalarType Ewald<ScalarType>::force_functor(
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
    template<class ScalarType>
    inline ScalarType Ewald<ScalarType>::pot_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = double(temp);
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const ScalarType interp = Internal::quadraticInterpolate<ScalarType>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return charges[i] * charges[j] * interp;
    }

    template<class ScalarType>
    void Ewald<ScalarType>::initIntegralLimit(ScalarType volume) {
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

    template<class ScalarType>
    void Ewald<ScalarType>::makeTables() {
        for (size_t i = 1; i < erfc_table.getLength(); ++i) {
            const auto x = PlainScalar((i - 1) * ErfcTableStep);
            erfc_table[i] = erfc(x) / x * integralLimit.getValue();
        }
        erfc_table[0] = erfc_table[1] = erfc_table[2]; // Smooth out divergent erfc(0) / 0 
        erfcStep = PlainScalar(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * PlainScalar(2));
    }

    template<class ScalarType>
    typename Ewald<ScalarType>::ComplexType
    Ewald<ScalarType>::kSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / 3;
        const size_t atom2 = dof2 / 3;
        const size_t direction1 = dof1 % 3U;
        const size_t direction2 = dof2 % 3U;
        const ScalarType charge1 = charges[atom1];
        const ScalarType charge2 = charges[atom2];

        const ScalarType factor1 = reciprocal(square(PlainScalar(2) * integralLimit));
        ComplexType kSpaceSum = 0;
        PeriodicCell<ScalarType, Dim>::forCellInRange(kSpaceSumRange, repCell.getLattice(),
            [this, atom1, atom2, direction1, direction2, charge1, charge2, waveQ, factor1, &pos, &kSpaceSum](Vector3D waveG) {
                ComplexType temp = 0;
                {
                    const Vector3D waveSum = waveQ + waveG;
                    const ScalarType squaredNorm = waveSum.squaredNorm();
                    const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                    if (isNotGammaPoint) {
                        const ScalarType factor2 = squaredNorm * exp(squaredNorm * factor1);
                        const ScalarType phase = waveSum * (pos.row(atom1).asVector() - pos.row(atom2).asVector());
                        const ComplexType elem = (waveSum[direction1] * waveSum[direction2]) / factor2 * ComplexType::fromPhase(phase);
                        temp = charge1 * charge2 * elem;
                    }
                }

                const bool isSameAtom = atom1 == atom2;
                if (isSameAtom) {
                    const ScalarType squaredNorm = waveG.squaredNorm();
                    const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                    if (isNotGammaPoint) {
                        const size_t numParticle = getNumParticle();
                        const Vector<ScalarType> dots = pos * waveG;
                        Vector<ScalarType> sin_vec(numParticle);
                        Vector<ScalarType> cos_vec(numParticle);
                        sincos(dots, sin_vec, cos_vec);
                        const ScalarType sum_cos = cos_vec * charges;
                        const ScalarType sum_sin = sin_vec * charges;
                        const ScalarType phase = waveG * pos.row(atom1).asVector();
                        const ComplexType elem = cos(phase) * sum_cos + sin(phase) * sum_sin;

                        const ScalarType factor2 = squaredNorm * exp(squaredNorm * factor1);
                        temp -= charge1 * elem * ((waveG[direction1] * waveG[direction2]) / factor2);
                    }
                }
                kSpaceSum += temp;
            });
        kSpaceSum *= PlainScalar(4 * M_PI) * inv_volume;
        return kSpaceSum;
    }

    template<class ScalarType>
    typename Ewald<ScalarType>::ComplexType
    Ewald<ScalarType>::rSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / 3;
        const size_t atom2 = dof2 / 3;
        const size_t direction1 = dof1 % 3U;
        const size_t direction2 = dof2 % 3U;
        const ScalarType charge1 = charges[atom1];
        const ScalarType charge2 = charges[atom2];

        ComplexType rSpaceSum = 0;
        PeriodicCell<ScalarType, Dim>::forCellInRange(rSpaceSumRange, lattice,
            [this, atom1, atom2, direction1, direction2, charge1, charge2, waveQ, &pos, &rSpaceSum](Vector3D delta) {
                const bool isSameDirection = direction1 == direction2;
                ComplexType temp = 0;
                {
                    const Vector3D x = pos.row(atom1).asVector() - pos.row(atom2).asVector() + delta;
                    const ScalarType norm = x.norm();
                    const bool isNotSelf = norm > ScalarType(std::numeric_limits<ScalarType>::min());
                    if (isNotSelf) {
                        const ScalarType phase = waveQ * delta;
                        temp = rSpaceForceConstImpl1(norm) * (x[direction1] * x[direction2]);
                        if (isSameDirection)
                            temp += rSpaceForceConstImpl2(norm);
                        temp *= ComplexType::fromPhase(phase) * (charge1 * charge2);
                    }
                }

                const bool isSameAtom = atom1 == atom2;
                if (isSameAtom) {
                    ComplexType temp1 = 0;
                    for (size_t i = 0; i < getNumParticle(); ++i) {
                        const Vector3D x = pos.row(atom1).asVector() - pos.row(i).asVector() + delta;
                        const ScalarType norm = x.norm();
                        const bool isNotSelf = norm > ScalarType(std::numeric_limits<ScalarType>::min());
                        if (isNotSelf) {
                            ComplexType temp2 = rSpaceForceConstImpl1(norm) * (x[direction1] * x[direction2]);
                            if (isSameDirection)
                                temp2 += rSpaceForceConstImpl2(norm);
                            temp1 += temp2 * charges[i];
                        }
                    }
                    temp -= temp1 * charge1;
                }
                rSpaceSum += temp;
            });
        return rSpaceSum;
    }

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::rSpaceForceConstImpl1(ScalarType r) const {
        const ScalarType x = integralLimit * r;
        const ScalarType x2 = square(x);
        const ScalarType term1 = ScalarType(3) * erfc(x);
        const ScalarType term2 = ScalarType(2 * M_2_SQRTPI) * x * (ScalarType(1.5) + x2) / exp(x2);
        return (term1 + term2) / (square(square(r)) * r);
    }

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::rSpaceForceConstImpl2(ScalarType r) const {
        const ScalarType x = integralLimit * r;
        const ScalarType term1 = erfc(x);
        const ScalarType term2 = ScalarType(M_2_SQRTPI) * x / exp(square(x));
        return -(term1 + term2) / (square(r) * r);
    }
}
