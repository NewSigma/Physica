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

#include "Physica/Core/Math/Calculus/ODE/SRK2.h"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Langevin.h"

namespace Physica::Core {
    /**
     * Note: Velocity rescaling algorithm[1] may give a negative probability at small number of particles condition.
     * 
     * Reference:
     * [1] G, Bussi, D. Donadio and M. Parrinello, J. Chem. Phys. 126, 014101 (2007).
     */
    template<class KineticModel>
    class DoubleThermo {
        using TraitsType = Traits<KineticModel>;
        using ScalarType = typename TraitsType::ScalarType;
        constexpr static unsigned int Dim = TraitsType::Dim;
        constexpr static size_t NumReplica = TraitsType::NumReplica;
        using MDCellType = MDCell<ScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, Dim, NumReplica>;
        using BufferType = typename RingPolymerType::BufferType;

        ScalarType temperatureT;
        ScalarType thermostatTime;
    public:
        DoubleThermo(ScalarType temperatureT_, ScalarType thermostatTime_);
        DoubleThermo(const DoubleThermo&) = default;
        DoubleThermo(DoubleThermo&&) noexcept = default;
        ~DoubleThermo() = default;
        /* Operators */
        DoubleThermo& operator=(DoubleThermo obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class RandomPoolType, class Executor>
        void step(RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const;
        void swap(DoubleThermo& __restrict obj) noexcept;
        /* Setters */
        void setTemperature(ScalarType temperatureT_) { temperatureT = temperatureT_; }
        void setThermostatTime(ScalarType time) { thermostatTime = time; }
    private:
        template<class RandomPoolType>
        ScalarType makeTranslationalFactor(const RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const;
    };

    template<class KineticModel>
    DoubleThermo<KineticModel>::DoubleThermo(ScalarType temperatureT_, ScalarType thermostatTime_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_) {}

    template<class KineticModel>
    template<class RandomPoolType, class Executor>
    void DoubleThermo<KineticModel>::step(
            RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const {
        const size_t dof = ringPolymer.getDOF();
        const ScalarType factor_translational = makeTranslationalFactor<RandomPoolType>(ringPolymer, deltaT, pool);
        if constexpr (NumReplica != 1) {
            const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
            const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
            auto future = Executor::parallel_for(
                [factor_translational, repBeta, omegaW, deltaT, &ringPolymer, &pool](unsigned int i) {
                    const size_t numReplica = ringPolymer.getNumReplica();
                    const auto& massVec = ringPolymer.getMassVec();

                    const auto mass = massVec[i / Dim];
                    const ScalarType factor = sqrt(repBeta * mass);
                    auto fft = FFT<ScalarType, 1>::makeEmptyFFT(numReplica);
                    BufferType buffer(2, ringPolymer.getKSpaceSize());

                    ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                    fft.getRSpace().random_normal(pool);
                    FFT<ScalarType, 1>::transform(ringPolymer.getFFT(), fft);
                    /* Translational mode */ {
                        buffer(0, 0) *= factor_translational;
                    }
                    for (size_t j = 1; j < buffer.getColumn(); ++j) {
                        const ScalarType phase = M_PI * j / numReplica;
                        const ScalarType viscosityY = sin(phase) * omegaW;
                        Langevin<ScalarType, Dim>::langevinImpl(
                                buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                    }
                    ringPolymer.toBeadRepr(i, ringPolymer.asMatrix(), buffer, fft);
                }, dof, Executor::getNumThread());
            Executor::auto_wait(future);
        }
        else {
            for (size_t i = 0; i < dof; ++i) {
                ringPolymer.asMatrix()(i, 0) *= factor_translational;
            }
        }
    }

    template<class KineticModel>
    void DoubleThermo<KineticModel>::swap(DoubleThermo& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
    }

    template<class KineticModel>
    template<class RandomPoolType>
    typename DoubleThermo<KineticModel>::ScalarType DoubleThermo<KineticModel>::makeTranslationalFactor(
            const RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const {
        using Integrator = SRK2<ScalarType, 1>;
        using VectorType = typename Integrator::VectorType;

        if (temperatureT.isZero())
            return 0;
        const size_t dof = ringPolymer.getDOF();
        [[maybe_unused]] ScalarType unused = 0;
        const ScalarType nowT = ringPolymer.template calcTemperature<KineticModel>();
        VectorType sol{nowT};
        Integrator::step(deltaT, unused, sol,
                            [this]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                return {(temperatureT - sol[0]) / thermostatTime};
                            },
                            [this, dof, &pool]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                const ScalarType rand = ScalarType::random_normal(pool);
                                return {sqrt((temperatureT * sol[0]) / (thermostatTime * dof)) * 2 * rand};
                            });
        if (!sol[0].isPositive()) [[unlikely]]
            throw std::invalid_argument("[Error]: Number of particle is too small that negative probability is encountered");
        return sqrt(temperatureT / sol[0]);
    }
}

namespace Physica {
    template<class KineticModel>
    class Traits<Core::DoubleThermo<KineticModel>> {
    public:
        constexpr static bool IsCentroidCoupled = false;
    };
}
