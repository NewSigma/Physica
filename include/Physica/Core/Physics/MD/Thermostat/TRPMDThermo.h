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

#include "Langevin.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Rossi M, Ceriotti M, Manolopoulos D E. How to remove the spurious resonances from ring polymer molecular dynamics[J]. J. Chem. Phys, 2014, 140(23):5106.
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class TRPMDThermo {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim, NumReplica>;
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
        template<class RandomGenerator>
        void step(RingPolymerType& ringPolymer, RandomGenerator& gen, ScalarType deltaT) const;
        void swap(TRPMDThermo& obj) noexcept;
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>::TRPMDThermo(ScalarType temperatureT_)
            : temperatureT(temperatureT_) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>&
    TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>::operator=(TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class RandomGenerator>
    void TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>::step(
            RingPolymerType& ringPolymer,
            RandomGenerator& gen,
            ScalarType deltaT) const {
        if constexpr (NumReplica == 1)
            return;
        const size_t dof = ringPolymer.getDOF();
        const size_t numReplica = ringPolymer.getNumReplica();
        const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
        const auto& massVec = ringPolymer.getMassVec();
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
                for (size_t j = 1; j < buffer.getColumn(); ++j) {
                    const ScalarType phase = M_PI * j / numReplica;
                    const ScalarType viscosityY = sin(phase) * omegaW;
                    Langevin<ScalarType, PosScalarType, Dim>::langevinImpl(
                            buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                }
                ringPolymer.toBeadRepr(i);
            }
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>::swap(TRPMDThermo<ScalarType, PosScalarType, Dim, NumReplica>& obj) noexcept {
        temperatureT.swap(obj.temperatureT);
    }
}
