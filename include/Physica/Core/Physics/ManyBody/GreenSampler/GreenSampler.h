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

#include "Physica/Core/Physics/ManyBody/DQMC.h"

namespace Physica {
    template<Scalar T>
    class GreenSampler {
        using This = GreenSampler<T>;
    protected:
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        const HubbardParams<T>& params;
    private:
        VectorND<Tv> signs;
        size_t cursor = 0;
    public:
        GreenSampler(const HubbardParams<T>& params, size_t numSample);
        GreenSampler(const This&) = delete;
        GreenSampler(This&&) noexcept = delete;
        ~GreenSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] Tv calcSign() const noexcept;

        void reset() { cursor = 0; }
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] auto getRepelU() const noexcept { return params.getRepelU(); }
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return params.getHoppingMatrix(); }

        [[nodiscard]] const auto& getSigns() const noexcept { return signs; }
        [[nodiscard]] size_t getNumSample() const noexcept { return signs.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    protected:
        void sample(T sign) noexcept;
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(const HubbardParams<T>& params, size_t numSample) : params(params), signs(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign() const noexcept -> Tv {
        return signs.mean();
    }

    template<Scalar T>
    void GreenSampler<T>::sample(T sign) noexcept {
        signs[cursor] = sign;
        cursor = (cursor + 1) % getNumSample();
    }
}
