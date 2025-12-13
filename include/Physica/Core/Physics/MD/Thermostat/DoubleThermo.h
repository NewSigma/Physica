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

#include "Physica/Core/Math/Calculus/ODE/SRK2.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Langevin.h"

namespace Physica {
    /**
     * Note: Velocity rescaling algorithm[1] may give a negative probability at small number of particles condition.
     * 
     * Reference:
     * [1] G, Bussi, D. Donadio and M. Parrinello, J. Chem. Phys. 126, 014101 (2007).
     */
    template<class KineticModel>
    class DoubleThermo {
        using This = DoubleThermo<KineticModel>;
        using TraitsType = Traits<KineticModel>;
        using T = TraitsType::ScalarType;
        constexpr static unsigned int Dim = TraitsType::Dim;
        constexpr static size_t NumReplica = TraitsType::NumReplica;
        using MDCellType = MDCell<T, Dim>;
        using MassVector = MDCellType::MassVector;
        using RingPolymerType = RingPolymer<T, Dim, NumReplica>;
        using BufferType = RingPolymerType::BufferType;

        T temperatureT;
        T thermostatTime;
    public:
        DoubleThermo() = default;
        DoubleThermo(T temperatureT_, T thermostatTime_);
        DoubleThermo(const This&) = default;
        DoubleThermo(This&&) noexcept = default;
        ~DoubleThermo() = default;
        /* Operators */
        This& operator=(DoubleThermo obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R, ExecutePolicy P>
        void step(RingPolymerType& ringPolymer, T deltaT) const;
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setTemperature(T temperatureT_) { temperatureT = temperatureT_; }
        void setThermostatTime(T time) { thermostatTime = time; }
    private:
        template<RNG R>
        T makeTranslationalFactor(const RingPolymerType& ringPolymer, T deltaT) const;
    };

    template<class KineticModel>
    DoubleThermo<KineticModel>::DoubleThermo(T temperatureT_, T thermostatTime_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_) {}

    template<class KineticModel>
    template<RNG R, ExecutePolicy P>
    void DoubleThermo<KineticModel>::step(
            RingPolymerType& ringPolymer, T deltaT) const {
        const size_t dof = ringPolymer.getDOF();
        const T factor_translational = makeTranslationalFactor<R>(ringPolymer, deltaT);
        if constexpr (NumReplica != 1) {
            const T repBeta = ringPolymer.calcRepBeta(temperatureT);
            const T omegaW = ringPolymer.calcOmegaW(temperatureT);
            parallel_for<P>(
                [factor_translational, repBeta, omegaW, deltaT, &ringPolymer](unsigned int i) {
                    const size_t numReplica = ringPolymer.getNumReplica();
                    const auto& massVec = ringPolymer.getMassVec();

                    const auto mass = massVec[i / Dim];
                    const T factor = sqrt(repBeta * mass);
                    auto fft = FFT<T, 1>::makeEmptyFFT(numReplica);
                    BufferType buffer(2, ringPolymer.getKSpaceSize());

                    ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                    fft.getRSpace().template random_normal<R>();
                    FFT<T, 1>::transform(ringPolymer.getFFT(), fft);
                    /* Translational mode */ {
                        buffer(0, 0) *= factor_translational;
                    }
                    for (size_t j = 1; j < buffer.getCol(); ++j) {
                        const T phase = M_PI * j / numReplica;
                        const T viscosityY = sin(phase) * omegaW;
                        Langevin<T, Dim>::langevinImpl(
                                buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                    }
                    ringPolymer.toBeadRepr(i, ringPolymer.asMatrix(), buffer, fft);
                }, dof, 0).wait();
        }
        else
            ringPolymer.asMatrix().col(0).head(dof) *= factor_translational;
    }

    template<class KineticModel>
    void DoubleThermo<KineticModel>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
    }

    template<class KineticModel>
    template<RNG R>
    auto DoubleThermo<KineticModel>::makeTranslationalFactor(
            const RingPolymerType& ringPolymer, T deltaT) const -> T {
        using Integrator = SRK2<T, 1>;
        using VectorType = Integrator::VectorType;

        if (temperatureT.isZero())
            return 0;
        const size_t dof = ringPolymer.getDOF();
        [[maybe_unused]] T unused = 0;
        const T nowT = ringPolymer.template calcTemperature<KineticModel>();
        VectorType sol{nowT};
        Integrator::step(deltaT, unused, sol,
                            [this]([[maybe_unused]] T x, VectorType sol) -> VectorType {
                                return {(temperatureT - sol[0]) / thermostatTime};
                            },
                            [this, dof]([[maybe_unused]] T x, VectorType sol) -> VectorType {
                                const T rand = T::template random_normal<R>();
                                return {sqrt((temperatureT * sol[0]) / (thermostatTime * dof)) * 2 * rand};
                            });
        if (!sol[0].isPositive()) [[unlikely]]
            throw std::invalid_argument("[Error]: Bad probability, maybe particle number is too small");
        return sqrt(temperatureT / sol[0]);
    }
}

namespace Physica {
    template<class KineticModel>
    class Traits<DoubleThermo<KineticModel>> {
    public:
        constexpr static bool IsCentroidCoupled = false;
    };
}
