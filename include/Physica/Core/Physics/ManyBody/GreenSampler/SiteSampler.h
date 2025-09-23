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
    public:
        enum Observable {
            Density,
            DoubleOccupy
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
        [[nodiscard]] T calcMeanWeighted() const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    };

    template<Scalar T>
    SiteSampler<T>::SiteSampler(size_t numSample) : Base(numSample), observes(numSample) {}

    template<Scalar T>
    template<Scalar U>
    void SiteSampler<T>::sample(const DQMC<U>& dqmc, Observable type) {
        T mean = 0;
        for (int i = 0; i < dqmc.getNumEqualGreen(); ++i) {
            const auto& greenU = dqmc.getGreenUs()[i];
            const auto& greenD = dqmc.getGreenDs()[i];

            T observe;
            switch (type) {
            case Density:
                observe = T(2) - greenU.diag().reals().mean() - greenD.diag().reals().mean();
                break;
            case DoubleOccupy:
                observe = (T(1) - greenU.diag().reals()) * (T(1) - greenD.diag().reals()) / T(dqmc.getNumSite());
                break;
            default:
                unreachable();
            }
            toNextMean(mean, i, observe);
        }
        observes[Base::getCursor()] = mean;
        Base::sample(dqmc.getLnPartitionZ(), dqmc.getSign());
    }

    template<Scalar T>
    T SiteSampler<T>::calcMean() const {
        return Base::calcMean(observes);
    }

    template<Scalar T>
    T SiteSampler<T>::calcMeanWeighted() const {
        return Base::calcMeanWeighted(observes);
    }

    template<Scalar T>
    void SiteSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        observes.swap(obj.observes);
    }
}
