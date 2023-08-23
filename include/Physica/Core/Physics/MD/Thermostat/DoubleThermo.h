/*
 * Copyright 2023 WeiBo He.
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
    template<class ScalarType, class PosScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class DoubleThermo {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim, NumReplica>;
        using BufferType = typename RingPolymerType::BufferType;

        ScalarType temperatureT;
        ScalarType thermostatTime;
    public:
        DoubleThermo(ScalarType temperatureT_, ScalarType thermostatTime_);
        DoubleThermo(const DoubleThermo&) = default;
        DoubleThermo(DoubleThermo&&) noexcept = default;
        ~DoubleThermo() = default;
        /* Operators */
        DoubleThermo& operator=(DoubleThermo obj) noexcept;
        /* Operations */
        template<class RandomGenerator, class Executor>
        void step(RingPolymerType& ringPolymer, ScalarType deltaT) const;
        void swap(DoubleThermo& obj) noexcept;
        /* Setters */
        void setThermostatTime(ScalarType time) { thermostatTime = time; }
    private:
        template<class RandomGenerator>
        ScalarType makeTranslationalFactor(const RingPolymerType& ringPolymer, ScalarType deltaT) const;
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>::DoubleThermo(ScalarType temperatureT_, ScalarType thermostatTime_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>&
    DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>::operator=(DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomGenerator, class Executor>
    void DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>::step(
            RingPolymerType& ringPolymer, ScalarType deltaT) const {
        const size_t dof = ringPolymer.getDOF();
        const ScalarType factor_translational = makeTranslationalFactor<RandomGenerator>(ringPolymer, deltaT);
        if constexpr (NumReplica != 1) {
            const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
            const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
            auto future = Executor::parallel_for([this, factor_translational, repBeta, omegaW, deltaT, &ringPolymer](unsigned int i) {
                const size_t numReplica = ringPolymer.getNumReplica();
                const auto& massVec = ringPolymer.getMassVec();

                const auto mass = massVec[i / Dim];
                const ScalarType factor = sqrt(repBeta * mass);
                auto fft = FFT<ScalarType, 1>::makeEmptyFFT(numReplica, 1);
                BufferType buffer(2, ringPolymer.getKSpaceSize());

                ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                fft.getRSpace().random_normal(RandomPool<RandomGenerator>::getGen());
                FFT<ScalarType, 1>::transform(ringPolymer.getFFT(), fft);
                /* Translational mode */ {
                    buffer(0, 0) *= factor_translational;
                }
                for (size_t j = 1; j < buffer.getColumn(); ++j) {
                    const ScalarType phase = M_PI * j / numReplica;
                    const ScalarType viscosityY = sin(phase) * omegaW;
                    Langevin<ScalarType, PosScalarType, Dim>::langevinImpl(
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

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>::swap(DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>& obj) noexcept {
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomGenerator>
    ScalarType DoubleThermo<ScalarType, PosScalarType, Dim, NumReplica>::makeTranslationalFactor(
            const RingPolymerType& ringPolymer, ScalarType deltaT) const {
        const size_t dof = ringPolymer.getDOF();
        using VectorType = Vector<ScalarType, 1>;
        ScalarType result = 1.0;
        [[maybe_unused]] ScalarType _ = 0;
        const ScalarType nowT = ringPolymer.calcTemperature();
        VectorType sol{nowT};
        SRK2<ScalarType, 1>::step(deltaT, _, sol,
                                    [this]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                        return {(temperatureT - sol[0]) / thermostatTime};
                                    },
                                    [this, dof]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                        std::normal_distribution dist{};
                                        RandomGenerator& gen = RandomPool<RandomGenerator>::getGen();
                                        return {sqrt((temperatureT * sol[0]) / (thermostatTime * dof)) * 2 * dist(gen)};
                                    });
        if (sol[0].isPositive()) [[likely]]
            result = sqrt(temperatureT / sol[0]);
        else
            throw std::invalid_argument("[Error]: Number of particle is too small that negative probability is encountered");
        return result;
    }
}
