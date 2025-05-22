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
        void sample(const DQMC<T>& dqmc, Observable type);
        [[nodiscard]] T calc() const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    };

    template<Scalar T>
    SiteSampler<T>::SiteSampler(size_t numSample) : Base(numSample), observes(numSample) {}

    template<Scalar T>
    void SiteSampler<T>::sample(const DQMC<T>& dqmc, Observable type) {
        T observe;
        switch (type) {
        case Density:
            observe = T(2) - dqmc.getGreenU().diag().mean() - dqmc.getGreenD().diag().mean();
            break;
        case DoubleOccupy:
            observe = (T(1) - dqmc.getGreenU().diag()) * (T(1) - dqmc.getGreenD().diag()) / T(dqmc.getNumSite());
            break;
        default:
            unreachable();
        }
        observes[Base::getCursor()] = observe;
        Base::sample(dqmc);
    }

    template<Scalar T>
    T SiteSampler<T>::calc() const {
        return Base::calc(observes);
    }

    template<Scalar T>
    void SiteSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        observes.swap(obj.observes);
    }
}
