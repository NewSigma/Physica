/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2023 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/Filter/SavitzkyGolay.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica {
    /**
     * Reference:
     * [1] J. Chem. Phys. 129, 074501 (2008); https://doi.org/10.1063/1.2968555
     */
    template<Scalar T>
    class IRAbsorb {
        using This = IRAbsorb<T>;
        using VectorType = VectorND<T>;
    public:
        VectorType dipoleCorr;
        FFT<T, 1> fft;
        T factor;
        SavitzkyGolay<T> filter;
    public:
        IRAbsorb(VectorType dipoleCorr_,
                 T temperatureT,
                 T deltaT,
                 T volume,
                 uint8_t filterRange,
                 size_t filterOrder);
        IRAbsorb(const IRAbsorb&) = default;
        IRAbsorb(IRAbsorb&&) noexcept = default;
        ~IRAbsorb() = default;
        /* Operators */
        IRAbsorb& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType makeWaveNum() const;
        [[nodiscard]] VectorType makeSpectrum();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDataLength() const noexcept { return dipoleCorr.getLength(); }
        [[nodiscard]] T getDeltaT() const noexcept { return filter.getDelta(); }
        [[nodiscard]] T getDeltaOmegaW() const noexcept { return fft.getKSpaceDelta(getDeltaT()); }
        [[nodiscard]] T getDeltaWaveNum() const noexcept { return getDeltaOmegaW() / T(2 * M_PI * PhyConst<AU>::speedOfLight); }
    };

    template<Scalar T>
    IRAbsorb<T>::IRAbsorb(VectorType dipoleCorr_,
                          T temperatureT,
                          T deltaT,
                          T volume,
                          uint8_t filterRange,
                          size_t filterOrder)
            : dipoleCorr(std::move(dipoleCorr_))
            , fft(dipoleCorr.getLength(), PlanFlag::Estimate)
            , factor(1 / (6 * PhyConst<AU>::vacuumDielectric * PhyConst<AU>::speedOfLight * double(volume) * PhyConst<AU>::boltzmannK * double(temperatureT)))
            , filter(filterRange, filterRange, filterOrder, deltaT) {}

    template<Scalar T>
    auto IRAbsorb<T>::makeWaveNum() const -> VectorType {
        const size_t size = fft.getKSpaceSize();
        return VectorType::linspace(0, getDeltaWaveNum() * (size - 1), size);
    }

    template<Scalar T>
    auto IRAbsorb<T>::makeSpectrum() -> VectorType {
        auto& rSpace = fft.getRSpace();
        rSpace = dipoleCorr;
        {
            for (size_t i = 0; i < rSpace.getLength(); ++i) {
                T hann = cos(T(M_PI * i / (rSpace.getLength() - 1))) * 0.5 + 0.5;
                rSpace[i] *= hann;
            }
            rSpace[0] *= T(0.5);
        }
        fft.transform();

        const auto& kSpace = fft.getKSpace();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const T deltaOmegaW = getDeltaOmegaW();
        VectorND<T> spectrum(kSpaceSize);
        for (size_t i = 0; i < kSpaceSize; ++i)
            spectrum[i] = factor * abs(kSpace[i].real()) * square(deltaOmegaW * i);
        filter.smooth_zero(spectrum);
        return spectrum.segment(filter.getLRange(), kSpaceSize + filter.getLRange());
    }

    template<Scalar T>
    void IRAbsorb<T>::swap(IRAbsorb& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        dipoleCorr.swap(obj.dipoleCorr);
        fft.swap(obj.fft);
        factor.swap(obj.factor);
        filter.swap(obj.filter);
    }
}
