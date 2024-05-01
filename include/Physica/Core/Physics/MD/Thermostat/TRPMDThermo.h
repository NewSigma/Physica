/*
 * Copyright 2023-2024 WeiBo He.
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

#include "Langevin.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic> class TRPMDThermo;

    namespace Internal {
        template<class ScalarType, unsigned int Dim, size_t NumReplica>
        class Traits<TRPMDThermo<ScalarType, Dim, NumReplica>> {
        public:
            constexpr static bool IsCentroidCoupled = false;
        };
    }
    /**
     * Reference:
     * [1] Rossi M, Ceriotti M, Manolopoulos D E. How to remove the spurious resonances from ring polymer molecular dynamics[J]. J. Chem. Phys, 2014, 140(23):5106.
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    class TRPMDThermo {
        using MDCellType = MDCell<ScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, Dim, NumReplica>;
        using BufferType = typename RingPolymerType::BufferType;

        ScalarType temperatureT;
    public:
        TRPMDThermo(ScalarType temperatureT_);
        TRPMDThermo(const TRPMDThermo&) = default;
        TRPMDThermo(TRPMDThermo&&) noexcept = default;
        ~TRPMDThermo() = default;
        /* Operators */
        TRPMDThermo& operator=(TRPMDThermo obj) noexcept;
        /* Operations */
        template<class RandomPoolType, class Executor>
        void step(RingPolymerType& ringPolymer, ScalarType deltaT, RandomPoolType& pool) const;
        void swap(TRPMDThermo& __restrict obj) noexcept;
        /* Setters */
        void setTemperature(ScalarType temperatureT_) { temperatureT = temperatureT_; }
    };

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    TRPMDThermo<ScalarType, Dim, NumReplica>::TRPMDThermo(ScalarType temperatureT_)
            : temperatureT(temperatureT_) {}

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    TRPMDThermo<ScalarType, Dim, NumReplica>&
    TRPMDThermo<ScalarType, Dim, NumReplica>::operator=(TRPMDThermo<ScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomPoolType, class Executor>
    void TRPMDThermo<ScalarType, Dim, NumReplica>::step(
            RingPolymerType& ringPolymer,
            ScalarType deltaT,
            RandomPoolType& pool) const {
        if constexpr (NumReplica == 1)
            return;
        const size_t dof = ringPolymer.getDOF();
        const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
        if constexpr (NumReplica != 1) {
            const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
            auto future = Executor::parallel_for(
                [this, repBeta, omegaW, deltaT, &ringPolymer, &pool](unsigned int i) {
                    const size_t numReplica = ringPolymer.getNumReplica();
                    const auto& massVec = ringPolymer.getMassVec();

                    const auto mass = massVec[i / Dim];
                    const ScalarType factor = sqrt(repBeta * mass);
                    auto fft = FFT<ScalarType, 1>::makeEmptyFFT(numReplica);
                    BufferType buffer(2, ringPolymer.getKSpaceSize());

                    ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                    fft.getRSpace().random_normal(pool);
                    FFT<ScalarType, 1>::transform(ringPolymer.getFFT(), fft);
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
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void TRPMDThermo<ScalarType, Dim, NumReplica>::swap(TRPMDThermo& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
    }
}
