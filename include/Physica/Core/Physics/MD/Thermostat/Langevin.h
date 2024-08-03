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

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;
    /**
     * PILE thermostat as introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
     */
    template<class ScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class Langevin {
        using MDCellType = MDCell<ScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, Dim, NumReplica>;
        using BufferType = typename RingPolymerType::BufferType;

        ScalarType temperatureT;
        ScalarType thermostatTime;
        bool removeDrift;
    public:
        Langevin(ScalarType temperatureT_, ScalarType thermostatTime_, bool removeDrift_);
        Langevin(const Langevin&) = default;
        Langevin(Langevin&&) noexcept = default;
        ~Langevin() = default;
        /* Operators */
        Langevin& operator=(Langevin obj) noexcept;
        /* Operations */
        template<class RandomPoolType, class Executor>
        void step(RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const;
        void swap(Langevin& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool isRemoveDriftEnabled() const noexcept { return removeDrift; }
        /* Setters */
        void setTemperature(ScalarType temperatureT_) { temperatureT = temperatureT_; }
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    Langevin<ScalarType, Dim, NumReplica>::Langevin(ScalarType temperatureT_, ScalarType thermostatTime_, bool removeDrift_)
            : temperatureT(temperatureT_)
            , thermostatTime(thermostatTime_)
            , removeDrift(removeDrift_) {}

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    Langevin<ScalarType, Dim, NumReplica>&
    Langevin<ScalarType, Dim, NumReplica>::operator=(Langevin<ScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomPoolType, class Executor>
    void Langevin<ScalarType, Dim, NumReplica>::step(
            RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const {
        const size_t dof = ringPolymer.getDOF();
        const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
        const ScalarType momentumViscosityY = Core::reciprocal(thermostatTime);
        const auto& massVec = ringPolymer.getMassVec();
        if constexpr (NumReplica != 1) {
            const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
            auto future = Executor::parallel_for(
                [deltaT, repBeta, omegaW, momentumViscosityY, &ringPolymer, &massVec, &pool](unsigned int i) {
                    const size_t numReplica = ringPolymer.getNumReplica();
                    const auto mass = massVec[i / Dim];
                    const ScalarType factor = sqrt(repBeta * mass);
                    auto fft = FFT<ScalarType, 1>::makeEmptyFFT(numReplica);
                    BufferType buffer(2, ringPolymer.getKSpaceSize());

                    ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                    fft.getRSpace().random_normal(pool);
                    FFT<ScalarType, 1>::transform(ringPolymer.getFFT(), fft);
                    /* Translational mode */ {
                        langevinImpl(buffer(0, 0), deltaT, momentumViscosityY, factor, fft.getKSpace()[0]);
                    }
                    for (size_t j = 1; j < buffer.getColumn(); ++j) {
                        const ScalarType phase = M_PI * j / numReplica;
                        const ScalarType viscosityY = sin(phase) * omegaW;
                        langevinImpl(buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                    }
                    ringPolymer.toBeadRepr(i, ringPolymer.asMatrix(), buffer, fft);
                }, dof, Executor::getNumThread());
            Executor::auto_wait(future);
        }
        else {
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                const ScalarType factor = sqrt(repBeta * mass);
                langevinImpl(ringPolymer.asMatrix()(i, 0), deltaT, momentumViscosityY, factor, ScalarType::random_normal(pool));
            }
        }

        if (removeDrift)
            ringPolymer.removeDrift();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class OtherScalar>
    void Langevin<ScalarType, Dim, NumReplica>::langevinImpl(
            OtherScalar& momentum,
            ScalarType deltaT,
            ScalarType viscosityY,
            ScalarType factor,
            OtherScalar random) {
        const ScalarType c1 = exp(-viscosityY * deltaT);
        const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
        momentum = c1 * momentum + factor * c2 * random;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void Langevin<ScalarType, Dim, NumReplica>::swap(Langevin& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
        thermostatTime.swap(obj.thermostatTime);
        std::swap(removeDrift, obj.removeDrift);
    }
}

namespace Physica {
    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    class Traits<Core::Langevin<ScalarType, Dim, NumReplica>> {
    public:
        constexpr static bool IsCentroidCoupled = true;
    };
}
