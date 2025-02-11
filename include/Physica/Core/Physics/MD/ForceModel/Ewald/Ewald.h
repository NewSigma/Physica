/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/PlainStruct.h"
#include "RSpaceEwald.h"

namespace Physica {
    /**
     * \tparam REwaldType: R(Space)EwaldType.
     * May be device side or host side, small cell version or large cell version, so it cannot be simplified to a single bool.
     * 
     * Reference:
     * [1] Martin, Richard M. Electronic structure: basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] Comput. Phys. Commun. 95(2-3), 73-92 (1996); https://doi.org/10.1016/0010-4655(96)00016-1
     */
    template<Scalar T, class REwaldType = RSpaceEwald<T, false>>
    class Ewald : private REwaldType {
        constexpr static bool IsDeviceREwald = is_device_obj<REwaldType>::value;
        using This = Ewald<T, REwaldType>;
        using Base = REwaldType;
        using Base::Dim;
        using typename Base::Tv;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using ComplexType = Complex<T>;
        using Index3D = GridBase::Index3D;
        using LatticeReturnType = std::conditional<IsDeviceREwald, LatticeMatrix, const LatticeMatrix&>::type;
        using HostChargeVector = std::conditional<IsDeviceREwald, VectorND<T>, PlainStruct<void>>::type;
    public:
        using typename Base::BornChargeArray;
    private:
        T selfE;
        T gammaPointE;
        [[no_unique_address]] HostChargeVector hostCharges;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice, VectorND<T> charges);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald obj) noexcept { swap(obj); return *this; }
        Ewald& operator=(Base base) noexcept;
        /* Operations */
        [[nodiscard]] T potentialV(const PositionMatrix& pos) const;

        template<class Executor>
        [[nodiscard]] VectorND<T> force(const PositionMatrix& pos);
        using Base::force_short;
        template<class Executor>
        [[nodiscard]] VectorND<T> force_long(const PositionMatrix& pos) const;

        [[nodiscard]] ComplexType forceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const;
        [[nodiscard]] ComplexType kSpaceForceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const;
        [[nodiscard]] ComplexType rSpaceForceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const;


        [[nodiscard]] LatticeMatrix virial(const PositionMatrix& pos);
        using Base::calcBornCharge;
        void swap(Ewald& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] inline LatticeReturnType getLattice() const noexcept;
        [[nodiscard]] inline LatticeReturnType getRepLattice() const noexcept;
        [[nodiscard]] inline const VectorND<T>& getCharges() const noexcept;
        using Base::getNumParticle;
        using Base::getVolume;
        using Base::getIntegralLimit;
        using Base::getKSpaceSumRange;
        [[nodiscard]] T getSelfE() const noexcept { return selfE; }
        [[nodiscard]] T getGammaPointE() const noexcept { return gammaPointE; }
        /* Setters */
        void setLattice(LatticeMatrix lattice);
    private:
        /* Operations */
        [[nodiscard]] ComplexType kSpaceForceConstImpl(const PositionMatrix& pos, const Vector3D<T>& waveQ, size_t dof1, size_t dof2) const;
        [[nodiscard]] inline ComplexType rSpaceForceConstImpl(const PositionMatrix& pos, const Vector3D<T>& waveQ, size_t dof1, size_t dof2) const;
        /* Friends */
        friend class ::Physica::Test;
    };

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::Ewald(LatticeMatrix lattice, VectorND<T> charges)
            : Base(std::move(lattice), std::move(charges)) {
        selfE = Base::calcSelfE();
        gammaPointE = Base::calcGammaPointE();
        if constexpr (IsDeviceREwald)
            hostCharges = std::move(charges);
    }

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>& Ewald<T, REwaldType>::operator=(Base base) noexcept {
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
    template<Scalar T, class REwaldType>
    T Ewald<T, REwaldType>::potentialV(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        const T rSpaceSum = Base::potentialV(pos);
        T kSpaceSum = 0;
        const T factor = reciprocal(square(Tv(2) * getIntegralLimit()));
        VectorND<T> dots(numParticle);
        VectorND<T> sin_vec(numParticle);
        VectorND<T> cos_vec(numParticle);
        PeriodicCell<T, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, factor, &pos, &kSpaceSum, &dots, &sin_vec, &cos_vec](Vector3D<T> delta) {
                const auto& charges = getCharges();
                const T squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const T sum_cos = cos_vec * charges;
                    const T sum_sin = sin_vec * charges;
                    kSpaceSum += (square(sum_cos) + square(sum_sin)) / (squaredNorm * exp(squaredNorm * factor));
                }
            });
        kSpaceSum *= Tv(4 * M_PI) * Base::getInvVolume();
        kSpaceSum += gammaPointE;
        return kSpaceSum + rSpaceSum - selfE;
    }

    template<Scalar T, class REwaldType>
    template<class Executor>
    VectorND<T> Ewald<T, REwaldType>::force(const PositionMatrix& pos) {
        VectorND<T> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            result = force_long<SeqExecutor>(pos);
        });
        const VectorND<T> rSpaceSum = Base::template force_short<Executor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }

    template<Scalar T, class REwaldType>
    template<class Executor>
    VectorND<T> Ewald<T, REwaldType>::force_long(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        VectorND<T> kSpaceSum(numParticle * Dim, 0);
        VectorND<T> dots(numParticle);
        VectorND<T> sin_vec(numParticle);
        VectorND<T> cos_vec(numParticle);
        const T factor1 = reciprocal(square(Tv(2) * getIntegralLimit()));
        PeriodicCell<T, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, numParticle, factor1, &dots, &sin_vec, &cos_vec, &pos, &kSpaceSum](Vector3D<T> delta) {
                const auto& charges = getCharges();
                const T squaredNorm = T(delta.squaredNorm());
                const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const T sum_cos = cos_vec * charges;
                    const T sum_sin = sin_vec * charges;
                    const T factor2 = reciprocal(squaredNorm * exp(squaredNorm * factor1));
                    for (size_t i = 0; i < numParticle; ++i) {
                        auto force_i = kSpaceSum.template segment<Dim>(i * Dim, (i + 1) * Dim);
                        const T temp = (sin_vec[i] * sum_cos - cos_vec[i] * sum_sin) * (factor2 * charges[i]);
                        force_i[0] += temp * delta[0];
                        force_i[1] += temp * delta[1];
                        force_i[2] += temp * delta[2];
                    }
                }
            });
        kSpaceSum *= T(8 * M_PI) * Base::getInvVolume(); // 8 Pi because time reversal symmetry
        return kSpaceSum;
    }
    /**
     * Reference:
     * [1] Rev. Mod. Phys. 73, 515; https://doi.org/10.1103/RevModPhys.73.515
     */
    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::ComplexType
    Ewald<T, REwaldType>::forceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const {
        const Vector3D<T> waveQ = getRepLattice().transpose() * qPoint;
        return kSpaceForceConstImpl(pos, waveQ, dof1, dof2) + rSpaceForceConstImpl(pos, waveQ, dof1, dof2);
    }

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::ComplexType
    Ewald<T, REwaldType>::kSpaceForceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const {
        const Vector3D<T> waveQ = getRepLattice().transpose() * qPoint;
        return kSpaceForceConstImpl(pos, waveQ, dof1, dof2);
    }

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::ComplexType
    Ewald<T, REwaldType>::rSpaceForceConst(const PositionMatrix& pos, Vector3D<T> qPoint, size_t dof1, size_t dof2) const {
        const Vector3D<T> waveQ = getRepLattice().transpose() * qPoint;
        return rSpaceForceConstImpl(pos, waveQ, dof1, dof2);
    }

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::LatticeMatrix Ewald<T, REwaldType>::virial(const PositionMatrix& pos) {
        const T factor = reciprocal(square(Tv(2) * getIntegralLimit()));
        const size_t numParticle = getNumParticle();
        VectorND<T> dots(numParticle);
        VectorND<T> sin_vec(numParticle);
        VectorND<T> cos_vec(numParticle);
        LatticeMatrix kSpaceVirial(Dim, Dim, 0);
        PeriodicCell<T, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            getKSpaceSumRange(), getRepLattice(), [this, &pos, &kSpaceVirial, &dots, &sin_vec, &cos_vec, factor](Vector3D<T> delta) {
                const auto& charges = getCharges();
                const T squaredNorm = T(delta.squaredNorm());
                const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                if (isNotGammaPoint) {
                    dots = pos * delta;
                    sincos(dots, sin_vec, cos_vec);
                    const T sum_cos = cos_vec * charges;
                    const T sum_sin = sin_vec * charges;
                    const T squaredStrucFactor = square(sum_cos) + square(sum_sin);

                    const T factor1 = squaredNorm * factor;
                    const T factor2 = squaredStrucFactor / (factor1 * exp(factor1));
                    const T factor3 = (factor1 + T(1)) * factor2 * (T(2) / squaredNorm);
                    LatticeMatrix temp = (factor3 * delta) * delta.transpose();
                    for (unsigned int i = 0; i < Dim; ++i)
                        temp(i, i) -= factor2;
                    kSpaceVirial += temp;
                }
            });
        kSpaceVirial *= T(-M_PI) * square(Base::getInvVolume() / getIntegralLimit());
        const T gammaPointP = gammaPointE * Base::getInvVolume();
        for (unsigned int i = 0; i < Dim; ++i)
            kSpaceVirial(i, i) += gammaPointP;
        const LatticeMatrix rSpaceVirial = Base::virial(pos);
        return kSpaceVirial + rSpaceVirial;
    }

    template<Scalar T, class REwaldType>
    void Ewald<T, REwaldType>::swap(Ewald& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        selfE.swap(obj.selfE);
        gammaPointE.swap(obj.gammaPointE);
        if constexpr (IsDeviceREwald)
            hostCharges.swap(obj.hostCharges);
    }

    template<Scalar T, class REwaldType>
    inline Ewald<T, REwaldType>::LatticeReturnType Ewald<T, REwaldType>::getLattice() const noexcept {
        if constexpr (IsDeviceREwald)
            return Base::getLattice().toHost();
        else
            return Base::getLattice();
    }

    template<Scalar T, class REwaldType>
    inline Ewald<T, REwaldType>::LatticeReturnType Ewald<T, REwaldType>::getRepLattice() const noexcept {
        if constexpr (IsDeviceREwald)
            return Base::getRepLattice().toHost();
        else
            return Base::getRepLattice();
    }

    template<Scalar T, class REwaldType>
    inline const VectorND<T>& Ewald<T, REwaldType>::getCharges() const noexcept {
        if constexpr (IsDeviceREwald)
            return hostCharges;
        else
            return Base::getCharges();
    }

    template<Scalar T, class REwaldType>
    void Ewald<T, REwaldType>::setLattice(LatticeMatrix lattice) {
        Base::setLattice(std::move(lattice));
        selfE = Base::calcSelfE();
        gammaPointE = Base::calcGammaPointE();
    }

    template<Scalar T, class REwaldType>
    Ewald<T, REwaldType>::ComplexType
    Ewald<T, REwaldType>::kSpaceForceConstImpl(const PositionMatrix& pos, const Vector3D<T>& waveQ, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / Dim;
        const size_t atom2 = dof2 / Dim;
        const size_t direction1 = dof1 % Dim;
        const size_t direction2 = dof2 % Dim;
        const auto& charges = getCharges();
        const T charge1 = charges[atom1];
        const T charge2 = charges[atom2];

        const T factor1 = reciprocal(square(Tv(2) * getIntegralLimit()));
        ComplexType kSpaceSum = 0;
        PeriodicCell<T, Dim>::forCellInRange(getKSpaceSumRange(), getRepLattice(),
            [this, atom1, atom2, direction1, direction2, charge1, charge2, waveQ, factor1, &pos, &kSpaceSum](Vector3D<T> waveG) {
                ComplexType temp = 0;
                {
                    const Vector3D<T> waveSum = waveQ + waveG;
                    const T squaredNorm = waveSum.squaredNorm();
                    const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                    if (isNotGammaPoint) {
                        const T factor2 = squaredNorm * exp(squaredNorm * factor1);
                        const T phase = waveSum * (pos.row(atom1) - pos.row(atom2));
                        const ComplexType elem = (waveSum[direction1] * waveSum[direction2]) / factor2 * ComplexType::fromPhase(phase);
                        temp = charge1 * charge2 * elem;
                    }
                }

                const bool isSameAtom = atom1 == atom2;
                if (isSameAtom) {
                    const auto& charges = getCharges();
                    const T squaredNorm = waveG.squaredNorm();
                    const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                    if (isNotGammaPoint) {
                        const size_t numParticle = getNumParticle();
                        const VectorND<T> dots = pos * waveG;
                        VectorND<T> sin_vec(numParticle);
                        VectorND<T> cos_vec(numParticle);
                        sincos(dots, sin_vec, cos_vec);
                        const T sum_cos = cos_vec * charges;
                        const T sum_sin = sin_vec * charges;
                        const T phase = waveG * pos.row(atom1);
                        const ComplexType elem = cos(phase) * sum_cos + sin(phase) * sum_sin;

                        const T factor2 = squaredNorm * exp(squaredNorm * factor1);
                        temp -= charge1 * elem * ((waveG[direction1] * waveG[direction2]) / factor2);
                    }
                }
                kSpaceSum += temp;
            });
        kSpaceSum *= Tv(4 * M_PI) * Base::getInvVolume();
        return kSpaceSum;
    }

    template<Scalar T, class REwaldType>
    inline Ewald<T, REwaldType>::ComplexType
    Ewald<T, REwaldType>::rSpaceForceConstImpl(const PositionMatrix& pos, const Vector3D<T>& waveQ, size_t dof1, size_t dof2) const {
        return Base::forceConst(pos, waveQ, dof1, dof2);
    }
}

namespace Physica {
    template<Scalar T1, class T2>
    class Traits<Ewald<T1, T2>> : public Traits<T2> {
    public:
        using ScalarType = T1;
        using REwaldType = T2;
    };
}
