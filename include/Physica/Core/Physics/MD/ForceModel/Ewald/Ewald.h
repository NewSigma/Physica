/*
 * Copyright 2021-2024 WeiBo He.
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

#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"
#include "Physica/Utils/CUDA/PlainStruct.h"
#include "RSpaceEwald.h"

namespace Physica::Core {
    template<class ScalarType, class REwaldType = RSpaceEwald<ScalarType, false>> class Ewald;

    namespace Internal {
        template<class T1, class T2>
        class Traits<Ewald<T1, T2>> : public Traits<T2> {
        public:
            using ScalarType = T1;
            using REwaldType = T2;
        };
    }
    /**
     * \tparam REwaldType: R(Space)EwaldType
     * 
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] Toukmaji A Y, Board J A. Ewald summation techniques in perspective: a survey[J]. Computer Physics Communications, 1996, 95(2-3):73-92.
     */
    template<class ScalarType, class REwaldType>
    class Ewald : private REwaldType {
        constexpr static bool IsDeviceREwald = Utils::is_device_obj<REwaldType>::value;
        using This = Ewald<ScalarType, REwaldType>;
        using Base = REwaldType;
        using Base::Dim;
        using typename Base::PlainScalar;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using ComplexType = ComplexScalar<ScalarType>;
        using CellListType = CellList<ScalarType>;
        using Index3D = typename CellListType::Index3D;
        using Vector3D = Vector<ScalarType, 3>;
        using LatticeReturnType = typename std::conditional<IsDeviceREwald, LatticeMatrix, const LatticeMatrix&>::type;
        using HostChargeVector = typename std::conditional<IsDeviceREwald, Vector<ScalarType>, PlainStruct<void>>::type;

        ScalarType selfE;
        ScalarType gammaPointE;
        HostChargeVector hostCharges;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice, Vector<ScalarType> charges);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald obj) noexcept { swap(obj); return *this; }
        Ewald& operator=(Base base) noexcept;
        /* Operations */
        [[nodiscard]] ScalarType potentialEnergy(const PositionMatrix& pos) const;

        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const PositionMatrix& pos);
        using Base::force_short;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const PositionMatrix& pos) const;

        [[nodiscard]] ComplexType forceConst(const PositionMatrix& pos, Vector3D qPoint, size_t dof1, size_t dof2) const;

        [[nodiscard]] LatticeMatrix virial(const PositionMatrix& pos);
        using Base::calcBornCharge;
        void swap(Ewald& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] inline LatticeReturnType getLattice() const noexcept;
        [[nodiscard]] inline LatticeReturnType getRepLattice() const noexcept;
        [[nodiscard]] inline const Vector<ScalarType>& getCharges() const noexcept;
        using Base::getNumParticle;
        using Base::getVolume;
        using Base::getIntegralLimit;
        using Base::getKSpaceSumRange;
        [[nodiscard]] ScalarType getSelfE() const noexcept { return selfE; }
        [[nodiscard]] ScalarType getGammaPointE() const noexcept { return gammaPointE; }
        /* Setters */
        void setLattice(LatticeMatrix lattice);
    private:
        /* Operations */
        [[nodiscard]] ComplexType kSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const;
        /* Friends */
        friend class ::Physica::Test;
    };

    template<class ScalarType, class REwaldType>
    Ewald<ScalarType, REwaldType>::Ewald(LatticeMatrix lattice, Vector<ScalarType> charges)
            : Base(std::move(lattice), std::move(charges)) {
        selfE = Base::calcSelfE();
        gammaPointE = Base::calcGammaPointE();
        if constexpr (IsDeviceREwald)
            hostCharges = std::move(charges);
    }

    template<class ScalarType, class REwaldType>
    Ewald<ScalarType, REwaldType>& Ewald<ScalarType, REwaldType>::operator=(Base base) noexcept {
        Base::swap(base);
        selfE = Base::calcSelfE();
        gammaPointE = Base::calcGammaPointE();
        if constexpr (IsDeviceREwald)
            Base::getCharges().toHost(hostCharges);
        return *this;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class ScalarType, class REwaldType>
    ScalarType Ewald<ScalarType, REwaldType>::potentialEnergy(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        const ScalarType rSpaceSum = Base::potentialEnergy(pos);
        ScalarType kSpaceSum = 0;
        const ScalarType factor = reciprocal(square(PlainScalar(2) * getIntegralLimit()));
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        PeriodicCell<ScalarType, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, numParticle, factor, &pos, &kSpaceSum, &dots, &sin_vec, &cos_vec](Vector3D delta) {
                const auto& charges = getCharges();
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const ScalarType sum_cos = cos_vec * charges;
                    const ScalarType sum_sin = sin_vec * charges;
                    kSpaceSum += (square(sum_cos) + square(sum_sin)) / (squaredNorm * exp(squaredNorm * factor));
                }
            });
        kSpaceSum *= PlainScalar(4 * M_PI) * Base::getInvVolume();
        kSpaceSum += gammaPointE;
        return kSpaceSum + rSpaceSum - selfE;
    }

    template<class ScalarType, class REwaldType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType, REwaldType>::force(const PositionMatrix& pos) {
        Vector<ScalarType> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            result = force_long<SequentialExecutor>(pos);
        });
        const Vector<ScalarType> rSpaceSum = Base::template force_short<Executor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }

    template<class ScalarType, class REwaldType>
    template<class Executor>
    Vector<ScalarType> Ewald<ScalarType, REwaldType>::force_long(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> kSpaceSum(numParticle * Dim, 0);
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        const ScalarType factor1 = reciprocal(square(PlainScalar(2) * getIntegralLimit()));
        PeriodicCell<ScalarType, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, numParticle, factor1, &dots, &sin_vec, &cos_vec, &pos, &kSpaceSum](Vector3D delta) {
                const auto& charges = getCharges();
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
        kSpaceSum *= ScalarType(8 * M_PI) * Base::getInvVolume(); // 8 Pi because time reversal symmetry
        return kSpaceSum;
    }
    /**
     * Reference:
     * [1] Rev. Mod. Phys. 73, 515; https://doi.org/10.1103/RevModPhys.73.515
     */
    template<class ScalarType, class REwaldType>
    typename Ewald<ScalarType, REwaldType>::ComplexType
    Ewald<ScalarType, REwaldType>::forceConst(const PositionMatrix& pos, Vector3D qPoint, size_t dof1, size_t dof2) const {
        const Vector3D waveQ = getRepLattice().transpose() * qPoint;
        const ComplexType kSpaceSum = kSpaceForceConst(pos, waveQ, dof1, dof2);
        const ComplexType rSpaceSum = Base::forceConst(pos, waveQ, dof1, dof2);
        return kSpaceSum + rSpaceSum;
    }

    template<class ScalarType, class REwaldType>
    typename Ewald<ScalarType, REwaldType>::LatticeMatrix Ewald<ScalarType, REwaldType>::virial(const PositionMatrix& pos) {
        const ScalarType factor = reciprocal(square(PlainScalar(2) * getIntegralLimit()));
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        LatticeMatrix kSpaceVirial(Dim, Dim, 0);
        PeriodicCell<ScalarType, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, &pos, &kSpaceVirial, &dots, &sin_vec, &cos_vec, factor](Vector3D delta) {
                const auto& charges = getCharges();
                const ScalarType squaredNorm = ScalarType(delta.squaredNorm());
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const ScalarType sum_cos = cos_vec * charges;
                    const ScalarType sum_sin = sin_vec * charges;
                    const ScalarType squaredStrucFactor = square(sum_cos) + square(sum_sin);

                    const ScalarType factor1 = squaredNorm * factor;
                    const ScalarType factor2 = squaredStrucFactor / (factor1 * exp(factor1));
                    const ScalarType factor3 = (factor1 + ScalarType(1)) * factor2 * (ScalarType(2) / squaredNorm);
                    LatticeMatrix temp = (factor3 * delta) * delta.transpose();
                    for (unsigned int i = 0; i < Dim; ++i)
                        temp(i, i) -= factor2;
                    kSpaceVirial += temp;
                }
            });
        kSpaceVirial *= ScalarType(-M_PI) * square(Base::getInvVolume() / getIntegralLimit());
        const ScalarType gammaPointP = gammaPointE * Base::getInvVolume();
        for (unsigned int i = 0; i < Dim; ++i)
            kSpaceVirial(i, i) += gammaPointP;
        const LatticeMatrix rSpaceVirial = Base::virial(pos);
        return kSpaceVirial + rSpaceVirial;
    }

    template<class ScalarType, class REwaldType>
    void Ewald<ScalarType, REwaldType>::swap(Ewald& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        selfE.swap(obj.selfE);
        gammaPointE.swap(obj.gammaPointE);
        if constexpr (IsDeviceREwald)
            hostCharges.swap(obj.hostCharges);
    }

    template<class ScalarType, class REwaldType>
    inline typename Ewald<ScalarType, REwaldType>::LatticeReturnType Ewald<ScalarType, REwaldType>::getLattice() const noexcept {
        if constexpr (IsDeviceREwald)
            return Base::getLattice().toHost();
        else
            return Base::getLattice();
    }

    template<class ScalarType, class REwaldType>
    inline typename Ewald<ScalarType, REwaldType>::LatticeReturnType Ewald<ScalarType, REwaldType>::getRepLattice() const noexcept {
        if constexpr (IsDeviceREwald)
            return Base::getRepLattice().toHost();
        else
            return Base::getRepLattice();
    }

    template<class ScalarType, class REwaldType>
    inline const Vector<ScalarType>& Ewald<ScalarType, REwaldType>::getCharges() const noexcept {
        if constexpr (IsDeviceREwald)
            return hostCharges;
        else
            return Base::getCharges();
    }

    template<class ScalarType, class REwaldType>
    void Ewald<ScalarType, REwaldType>::setLattice(LatticeMatrix lattice) {
        Base::setLattice(std::move(lattice));
        selfE = Base::calcSelfE();
        gammaPointE = Base::calcGammaPointE();
    }

    template<class ScalarType, class REwaldType>
    typename Ewald<ScalarType, REwaldType>::ComplexType
    Ewald<ScalarType, REwaldType>::kSpaceForceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / 3;
        const size_t atom2 = dof2 / 3;
        const size_t direction1 = dof1 % 3U;
        const size_t direction2 = dof2 % 3U;
        const auto& charges = getCharges();
        const ScalarType charge1 = charges[atom1];
        const ScalarType charge2 = charges[atom2];

        const ScalarType factor1 = reciprocal(square(PlainScalar(2) * getIntegralLimit()));
        ComplexType kSpaceSum = 0;
        PeriodicCell<ScalarType, Dim>::forCellInRange(getKSpaceSumRange(), getRepLattice(),
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
                    const auto& charges = getCharges();
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
        kSpaceSum *= PlainScalar(4 * M_PI) * Base::getInvVolume();
        return kSpaceSum;
    }
}
