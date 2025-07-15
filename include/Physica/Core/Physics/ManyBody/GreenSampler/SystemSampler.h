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
        using MatrixND = DenseMatrix<T>;
    public:
        enum Observable {
            AFM,  // Antiferromagnetic Structure Factor
            CDW,  // Charge Density Wave
            DOW   // Double Occupancy Wave
        };
    private:
        Array<MatrixND> observes;
        FFT<T, 2> fft;
    public:
        SystemSampler(const LatticeModel<2>& lattice, size_t numSample);
        SystemSampler(const This&) = default;
        SystemSampler(This&&) noexcept = default;
        ~SystemSampler() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        void sample(const DQMC<T>& dqmc, Observable type);
        [[nodiscard]] MatrixND calc() const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumSample;
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
        [[nodiscard]] int getNumSiteX() const noexcept { return fft.getRSpaceSize()[0]; }
        [[nodiscard]] int getNumSiteY() const noexcept { return fft.getRSpaceSize()[1]; }
        [[nodiscard]] int getNumSite() const noexcept { return getNumSiteX() * getNumSiteY(); }
    };

    template<Scalar T>
    SystemSampler<T>::SystemSampler(const LatticeModel<2>& lattice, size_t numSample)
            : Base(numSample), observes(numSample), fft(lattice.getSuperSize(), PlanFlag::Estimate) {}

    template<Scalar T>
    void SystemSampler<T>::sample(const DQMC<T>& dqmc, Observable type) {
        const int numSite = getNumSite();
        assert(numSite == dqmc.getNumSite() && "[Error]: Inconsistent site numbers");
        auto flatten = fft.getRSpace().flatten();
        for (int i = 0; i < dqmc.getNumEqualGreen(); ++i) {
            const auto& greenU = dqmc.getGreenUs()[i];
            const auto& greenD = dqmc.getGreenDs()[i];
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
            toNextMean(observes[Base::getCursor()], i, fft.getKSpace().squaredNorms() * reciprocal(T(numSite)));
        }
        Base::sample(dqmc);
    }

    template<Scalar T>
    auto SystemSampler<T>::calc() const -> MatrixND {
        const int kX = getNumSiteX();
        const int kY = FFT<T, 1>::rSizeToKSize(getNumSiteY());
        MatrixND result(kX, kY);
        VectorND<T> buffer(getNumSample());
        for (int x = 0; x < kX; ++x) {
            for (int y = 0; y < kY; ++y) {
                for (size_t i = 0; i < observes.getLength(); ++i)
                    buffer[i] = observes[i](x, y);
                result(x, y) = Base::calc(buffer);
            }
        }
        return result;
    }

    template<Scalar T>
    void SystemSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        observes.swap(obj.observes);
    }
}
