/*
 * Copyright 2025 Weibo He.
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

#include "GreenSampler.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica {
    /**
     * Defination of AFM, refer to [1]
     *
     * Reference:
     * [1] Phys. Rev. Lett. 62, 591; https://doi.org/10.1103/PhysRevLett.62.591
     */
    template<Scalar T>
    class SystemSampler : public GreenSampler<T> {
        using This = SystemSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = DQMC<T>::GreenPair;

        using typename Base::Tv;
    public:
        enum Observable {
            AFM, // Antiferromagnetic Structure Factor
            CDW, // Charge Density Wave
            DOW // Double Occupancy Wave
        };
    private:
        Array<MatrixND<T>> observes;
        FFT<T, 2> fft;
    public:
        SystemSampler(const DQMC<T>& dqmc, const LatticeModel<2>& lattice, size_t numSample);
        SystemSampler(const This&) = delete;
        SystemSampler(This&&) noexcept = delete;
        ~SystemSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const GreenPair& greens, Observable type);

        [[nodiscard]] MatrixND<T> calcMean() const;
        /* Getters */
        using Base::getNumSample;
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
        [[nodiscard]] int getNumSiteX() const noexcept { return fft.getRSpaceSize()[0]; }
        [[nodiscard]] int getNumSiteY() const noexcept { return fft.getRSpaceSize()[1]; }
        [[nodiscard]] int getNumSite() const noexcept { return getNumSiteX() * getNumSiteY(); }
    private:
        MatrixND<T> calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable type) noexcept;
    };

    template<Scalar T>
    SystemSampler<T>::SystemSampler(const DQMC<T>& dqmc, const LatticeModel<2>& lattice, size_t numSample)
            : Base(dqmc, numSample)
            , observes(numSample)
            , fft(lattice.getSuperSize(), PlanFlag::Estimate) {
        for (auto& elem : observes)
            elem.resize(fft.getKSpace());
    }

    template<Scalar T>
    void SystemSampler<T>::sample(const GreenPair& greens, Observable type) {
        observes[Base::getCursor()] = calcObservable(greens[0], greens[1], type) * Base::dqmc.getSign();
        Base::sample();
    }

    template<Scalar T>
    auto SystemSampler<T>::calcMean() const -> MatrixND<T> {
        const int kX = getNumSiteX();
        const int kY = FFT<T, 1>::rSizeToKSize(getNumSiteY());
        const Tv sign = Base::calcSign();
        MatrixND<T> result(kX, kY);
        VectorND<T> buffer(getNumSample());
        for (int x = 0; x < kX; ++x) {
            for (int y = 0; y < kY; ++y) {
                for (size_t i = 0; i < observes.getLength(); ++i)
                    buffer[i] = observes[i](x, y);
                result(x, y) = buffer.mean() / sign;
            }
        }
        return result;
    }

    template<Scalar T>
    auto SystemSampler<T>::calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable type) noexcept -> MatrixND<T> {
        int numSite = Base::getNumSite();
        auto flatten = fft.getRSpace().flatten();
        switch (type) {
        case AFM:
            for (int i = 0; i < numSite; ++i)
                flatten[i] = greenU.diag()[i] - greenD.diag()[i];
            break;
        case CDW:
            for (int i = 0; i < numSite; ++i)
                flatten[i] = T(2) - greenU.diag()[i] - greenD.diag()[i];
            break;
        case DOW:
            for (int i = 0; i < numSite; ++i)
                flatten[i] = (T(1) - greenU.diag()[i]) * (T(1) - greenD.diag()[i]);
            break;
        default:
            unreachable();
        }
        fft.transform();
        return fft.getKSpace().squaredNorms() * reciprocal(T(numSite));
    }
}
