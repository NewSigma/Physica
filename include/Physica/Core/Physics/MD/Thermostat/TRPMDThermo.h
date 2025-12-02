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

#include "Langevin.h"

namespace Physica {
    /**
     * Reference:
     * [1] J. Chem. Phys. 140, 234116 (2014); https://doi.org/10.1063/1.4883861
     */
    template<Scalar T, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class TRPMDThermo {
        using This = TRPMDThermo<T, Dim, NumReplica>;
        using MDCellType = MDCell<T, Dim>;
        using MassVector = MDCellType::MassVector;
        using RingPolymerType = RingPolymer<T, Dim, NumReplica>;
        using BufferType = RingPolymerType::BufferType;

        T temperatureT;
    public:
        TRPMDThermo() = default;
        TRPMDThermo(T temperatureT_);
        TRPMDThermo(const This&) = default;
        TRPMDThermo(This&&) noexcept = default;
        ~TRPMDThermo() = default;
        /* Operators */
        This& operator=(This obj) noexcept{ swap(obj); return *this; }
        /* Operations */
        template<RNG R, ExecutePolicy P>
        void step(RingPolymerType& ringPolymer, T deltaT) const;
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setTemperature(T temperatureT_) noexcept { temperatureT = temperatureT_; }
    };

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    TRPMDThermo<T, Dim, NumReplica>::TRPMDThermo(T temperatureT_)
            : temperatureT(temperatureT_) {}

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<RNG R, ExecutePolicy P>
    void TRPMDThermo<T, Dim, NumReplica>::step(RingPolymerType& ringPolymer, T deltaT) const {
        if constexpr (NumReplica == 1)
            return;
        const size_t dof = ringPolymer.getDOF();
        const T repBeta = ringPolymer.calcRepBeta(temperatureT);
        if constexpr (NumReplica != 1) {
            const T omegaW = ringPolymer.calcOmegaW(temperatureT);
            parallel_for<P>([repBeta, omegaW, deltaT, &ringPolymer](unsigned int i) {
                const size_t numReplica = ringPolymer.getNumReplica();
                const auto& massVec = ringPolymer.getMassVec();

                const auto mass = massVec[i / Dim];
                const T factor = sqrt(repBeta * mass);
                auto fft = FFT<T, 1>::makeEmptyFFT(numReplica);
                BufferType buffer(2, ringPolymer.getKSpaceSize());

                ringPolymer.toNormalRepr(i, ringPolymer.asMatrix(), buffer, fft);
                fft.getRSpace().template random_normal<R>();
                FFT<T, 1>::transform(ringPolymer.getFFT(), fft);
                for (size_t j = 1; j < buffer.getCol(); ++j) {
                    const T phase = M_PI * j / numReplica;
                    const T viscosityY = sin(phase) * omegaW;
                    Langevin<T, Dim>::langevinImpl(buffer(0, j), deltaT, viscosityY, factor, fft.getKSpace()[j]);
                }
                ringPolymer.toBeadRepr(i, ringPolymer.asMatrix(), buffer, fft);
            }, dof, 0).wait();
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void TRPMDThermo<T, Dim, NumReplica>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        temperatureT.swap(obj.temperatureT);
    }
}

namespace Physica {
    template<Scalar T, unsigned int Dim, size_t NumReplica>
    class Traits<TRPMDThermo<T, Dim, NumReplica>> {
    public:
        constexpr static bool IsCentroidCoupled = false;
    };
}
