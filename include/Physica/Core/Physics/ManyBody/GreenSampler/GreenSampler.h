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

        const DQMC<T>& dqmc;
    private:
        VectorND<Tr> lnAbsDets;
        VectorND<Tv> signs;
        size_t cursor = 0;
    public:
        GreenSampler(const DQMC<T>& dqmc_, size_t numSample);
        GreenSampler(const This&) = delete;
        GreenSampler(This&&) noexcept = delete;
        ~GreenSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T lnPartitionZ() const;
        [[nodiscard]] Tv calcSign() const noexcept;

        void reset() { cursor = 0; }
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return dqmc.getNumSite(); }
        [[nodiscard]] auto getRepelU() const noexcept { return dqmc.getParams().getRepelU(); }
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return dqmc.getParams().getHoppingMatrix(); }

        [[nodiscard]] const auto& getLnAbsDets() const noexcept { return lnAbsDets; }
        [[nodiscard]] const auto& getSigns() const noexcept { return signs; }
        [[nodiscard]] size_t getNumSample() const noexcept { return signs.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    protected:
        void sample() noexcept;
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(const DQMC<T>& dqmc_, size_t numSample) : dqmc(dqmc_), lnAbsDets(numSample), signs(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    T GreenSampler<T>::lnPartitionZ() const {
        const Tr factor = lnAbsDets.max();
        return ln(exp(lnAbsDets - factor) * getSigns()) + factor;
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign() const noexcept -> Tv {
        return signs.mean();
    }

    template<Scalar T>
    void GreenSampler<T>::sample() noexcept {
        lnAbsDets[cursor] = dqmc.getLnAbsDet();
        signs[cursor] = dqmc.getSign();
        cursor = (cursor + 1) % getNumSample();
    }
}
