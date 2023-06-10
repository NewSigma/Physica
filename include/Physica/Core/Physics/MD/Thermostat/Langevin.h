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

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;
    /**
     * Reference:
     * [1] M. Ceriotti, M. Parrinello, T. E. Markland and D. E. Manolopoulos, J. Chem. Phys. 133, 124104 (2010).
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class Langevin {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim, NumReplica>;
        using BufferType = typename RingPolymerType::BufferType;

        ScalarType temperatureT;
        ScalarType thermostatTime;
    public:
        Langevin(ScalarType temperatureT_, ScalarType thermostatTime_);
        Langevin(const Langevin&) = default;
        Langevin(Langevin&&) noexcept = default;
        ~Langevin() = default;
        /* Operators */
        Langevin& operator=(Langevin obj) noexcept;
        /* Operations */
        template<class RandomGenerator>
        void step(RingPolymerType& ringPolymer, RandomGenerator& gen, ScalarType deltaT) const;
        void swap(Langevin& obj) noexcept;
        /* Setters */
        void setThermostatTime(ScalarType time) { thermostatTime = time; }
        /* Static members */
        template<class OtherScalar>
        static inline void langevinImpl(
            OtherScalar& momentum,
            ScalarType deltaT,
            ScalarType viscosityY,
            ScalarType factor,
            OtherScalar random);
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    Langevin<ScalarType, PosScalarType, Dim, NumReplica>::Langevin(ScalarType temperatureT_, ScalarType thermostatTime_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    Langevin<ScalarType, PosScalarType, Dim, NumReplica>&
    Langevin<ScalarType, PosScalarType, Dim, NumReplica>::operator=(Langevin<ScalarType, PosScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomGenerator>
    void Langevin<ScalarType, PosScalarType, Dim, NumReplica>::step(
            RingPolymerType& ringPolymer,
            RandomGenerator& gen,
            ScalarType deltaT) const {
        const size_t dof = ringPolymer.getDOF();
        const size_t numReplica = ringPolymer.getNumReplica();
        const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
        const ScalarType momentumViscosityY = Core::reciprocal(thermostatTime);
        const auto& massVec = ringPolymer.getMassVec();
        std::normal_distribution<> dist{};
        if constexpr (NumReplica != 1) {
            const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
            FFT<ScalarType, 1> fft(numReplica, 1);
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                const ScalarType factor = sqrt(repBeta * mass);

                auto& buffer = ringPolymer.getBuffer();
                ringPolymer.toNormalRepr(i);

                fft.getRSpace().random_normal(gen);
                fft.transform();
                /* Translational mode */ {
                    langevinImpl(buffer(0, 0), deltaT, momentumViscosityY, factor, fft.getKSpace()[0]);
                }
                for (size_t j = 1; j < buffer.getColumn(); ++j) {
                    const ScalarType phase = M_PI * j / numReplica;
                    const ScalarType viscosityY = sin(phase) * omegaW;
                    langevinImpl(buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                }
                ringPolymer.toBeadRepr(i);
            }
        }
        else {
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                const ScalarType factor = sqrt(repBeta * mass);
                langevinImpl(ringPolymer.asMatrix()(i, 0), deltaT, momentumViscosityY, factor, ScalarType(dist(gen)));
            }
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class OtherScalar>
    void Langevin<ScalarType, PosScalarType, Dim, NumReplica>::langevinImpl(
            OtherScalar& momentum,
            ScalarType deltaT,
            ScalarType viscosityY,
            ScalarType factor,
            OtherScalar random) {
        const ScalarType c1 = exp(-viscosityY * deltaT);
        const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
        momentum = c1 * momentum + factor * c2 * random;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void Langevin<ScalarType, PosScalarType, Dim, NumReplica>::swap(Langevin<ScalarType, PosScalarType, Dim, NumReplica>& obj) noexcept {
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
    }
}
