/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica {
    template<Scalar T, unsigned int Dim, size_t NumReplica> class RingPolymer;
    /**
     * PILE thermostat as introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
     */
    template<Scalar T, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class Langevin {
        using This = Langevin<T, Dim, NumReplica>;
        using MDCellType = MDCell<T, Dim>;
        using MassVector = MDCellType::MassVector;
        using RingPolymerType = RingPolymer<T, Dim, NumReplica>;
        using BufferType = RingPolymerType::BufferType;

        T temperatureT;
        T thermostatTime;
        bool removeDrift;
    public:
        Langevin() = default;
        Langevin(T temperatureT_, T thermostatTime_, bool removeDrift_);
        Langevin(const This&) = default;
        Langevin(This&&) noexcept = default;
        ~Langevin() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R, ExecutePolicy P>
        void step(RingPolymerType& ringPolymer, T deltaT) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool isRemoveDriftEnabled() const noexcept { return removeDrift; }
        /* Setters */
        void setTemperature(T temperatureT_) { temperatureT = temperatureT_; }
        void setThermostatTime(T time) { thermostatTime = time; }
        /* Static members */
        template<Scalar U>
        static inline void langevinImpl(
            U& momentum,
            T deltaT,
            T viscosityY,
            T factor,
            U random);
    };

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    Langevin<T, Dim, NumReplica>::Langevin(T temperatureT_, T thermostatTime_, bool removeDrift_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_)
            , removeDrift(removeDrift_) {}

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<RNG R, ExecutePolicy P>
    void Langevin<T, Dim, NumReplica>::step(RingPolymerType& ringPolymer, T deltaT) const {
        const size_t dof = ringPolymer.getDOF();
        const T repBeta = ringPolymer.calcRepBeta(temperatureT);
        const T momentumViscosityY = reciprocal(thermostatTime);
        const auto& massVec = ringPolymer.getMassVec();
        if constexpr (NumReplica != 1) {
            const T omegaW = ringPolymer.calcOmegaW(temperatureT);
            parallel_for<P>(
                [deltaT, repBeta, omegaW, momentumViscosityY, &ringPolymer, &massVec](unsigned int i) {
                    const size_t numReplica = ringPolymer.getNumReplica();
                    const auto mass = massVec[i / Dim];
                    const T factor = sqrt(repBeta * mass);
                    auto fft = FFT<T, 1>::makeEmptyFFT(numReplica);
                    BufferType buffer(2, ringPolymer.getKSpaceSize());

                    ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                    fft.getRSpace().template random_normal<R>();
                    FFT<T, 1>::transform(ringPolymer.getFFT(), fft);
                    /* Translational mode */ {
                        langevinImpl(buffer(0, 0), deltaT, momentumViscosityY, factor, fft.getKSpace()[0]);
                    }
                    for (size_t j = 1; j < buffer.getCol(); ++j) {
                        const T phase = M_PI * j / numReplica;
                        const T viscosityY = sin(phase) * omegaW;
                        langevinImpl(buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                    }
                    ringPolymer.toBeadRepr(i, ringPolymer.asMatrix(), buffer, fft);
                }, dof, 0).wait();
        }
        else {
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                const T factor = sqrt(repBeta * mass);
                langevinImpl(ringPolymer.asMatrix()(i, 0), deltaT, momentumViscosityY, factor, T::template random_normal<R>());
            }
        }

        if (removeDrift)
            ringPolymer.removeDrift();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<Scalar U>
    void Langevin<T, Dim, NumReplica>::langevinImpl(
            U& momentum,
            T deltaT,
            T viscosityY,
            T factor,
            U random) {
        const T c1 = exp(-viscosityY * deltaT);
        const T c2 = sqrt(T(1) - square(c1));
        momentum = c1 * momentum + factor * c2 * random;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void Langevin<T, Dim, NumReplica>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
        std::swap(removeDrift, obj.removeDrift);
    }
}

namespace Physica {
    template<Scalar T, unsigned int Dim, size_t NumReplica>
    class Traits<Langevin<T, Dim, NumReplica>> {
    public:
        constexpr static bool IsCentroidCoupled = true;
    };
}
