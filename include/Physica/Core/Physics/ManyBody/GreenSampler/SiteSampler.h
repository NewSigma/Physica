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

namespace Physica {
    template<Scalar T>
    class SiteSampler : public GreenSampler<T> {
        using This = SiteSampler<T>;
        using Base = GreenSampler<T>;
        using typename Base::Trv;
    public:
        enum Observable : char {
            Density,
            DoubleOccupy,
            Kinetic,
            Potential,
            Internal
        };
    private:
        VectorND<T> observes;
    public:
        SiteSampler(size_t numSample);
        SiteSampler(const This&) = default;
        SiteSampler(This&&) noexcept = default;
        ~SiteSampler() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<Scalar U>
        void sample(const DQMC<U>& dqmc, Observable type);
        [[nodiscard]] T calcMean() const;
        [[nodiscard]] T calcMean(const VectorND<T>& lnWeights) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    private:
        template<Scalar U>
        static T calcObservable(const DQMC<U>& dqmc, int split, Observable type) noexcept;
    };

    template<Scalar T>
    SiteSampler<T>::SiteSampler(size_t numSample)
            : Base(numSample), observes(numSample) {}

    template<Scalar T>
    template<Scalar U>
    void SiteSampler<T>::sample(const DQMC<U>& dqmc, Observable type) {
        T mean = 0;
        for (int i = 0; i < dqmc.getNumEqualGreen(); ++i)
            mean.toNextMean(i, calcObservable(dqmc, i, type));
        observes[Base::getCursor()] = mean;
        Base::sample(dqmc.getLnAbsDet(), dqmc.getSign());
    }

    template<Scalar T>
    T SiteSampler<T>::calcMean() const {
        return Base::calcMean(observes);
    }

    template<Scalar T>
    T SiteSampler<T>::calcMean(const VectorND<T>& lnWeights) const {
        return Base::calcMean(observes, lnWeights);
    }

    template<Scalar T>
    void SiteSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        observes.swap(obj.observes);
    }

    template<Scalar T>
    template<Scalar U>
    T SiteSampler<T>::calcObservable(const DQMC<U>& dqmc, int split, Observable type) noexcept {
        const auto& greenU = dqmc.getGreenUs()[split];
        const auto& greenD = dqmc.getGreenDs()[split];
        switch (type) {
        case Density:
            return T(2) - greenU.diag().reals().mean() - greenD.diag().reals().mean();
        case DoubleOccupy:
            return (T(1) - greenU.diag().reals()) * (T(1) - greenD.diag().reals()) / T(dqmc.getNumSite());
        case Kinetic: {
            const auto& hoppingT = dqmc.getParams().getHoppingMatrix();
            return -hadamard(hoppingT, greenU + greenD).sum().real() / T(dqmc.getNumSite());
        }
        case Potential:
            return calcObservable(dqmc, split, DoubleOccupy) * dqmc.getParams().getRepelU();
        case Internal:
            return calcObservable(dqmc, split, Kinetic) + calcObservable(dqmc, split, Potential);
        default:
            unreachable();
        }
    }
}
