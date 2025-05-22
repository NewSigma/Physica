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
        using Tf = Diff<T, DiffMode::Forward, 1>;
    private:
        VectorND<Tf> samples;
        size_t cursor = 0;
    public:
        GreenSampler(size_t numSample);
        GreenSampler(const This&) = default;
        GreenSampler(This&&) noexcept = default;
        ~GreenSampler() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        void sample(const DQMC<T>& dqmc);
        [[nodiscard]] T calc(const VectorND<T>& observes) const;
        [[nodiscard]] T calcSign() const;
        [[nodiscard]] T lnPartitionZ() const;

        void reset() { cursor = 0; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getLnPartitionZs() const noexcept { return samples.values(); }
        [[nodiscard]] const auto& getSigns() const noexcept { return samples.grads(); }
        [[nodiscard]] size_t getNumSample() const noexcept { return samples.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(size_t numSample) : samples(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    void GreenSampler<T>::sample(const DQMC<T>& dqmc) {
        samples.values()[cursor] = dqmc.getLnPartitionZ();
        samples.grads()[cursor] = dqmc.getSign();
        cursor = (cursor + 1) % getNumSample();

        bool full = cursor == 0;
        if (full)
            samples -= getLnPartitionZs().max();
    }

    template<Scalar T>
    T GreenSampler<T>::calc(const VectorND<T>& observes) const {
        VectorND<Tf> buffer(getLnPartitionZs(), observes);
        const Tf sum = exp(buffer) * getSigns();
        return sum.grad() / sum.value(); // Avoid ln, its grad is well defined but value may give NAN
    }

    template<Scalar T>
    T GreenSampler<T>::calcSign() const {
        return samples.lnSumExp().grad();
    }

    template<Scalar T>
    T GreenSampler<T>::lnPartitionZ() const {
        return ln(exp(getLnPartitionZs()) * getSigns());
    }

    template<Scalar T>
    void GreenSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        std::swap(cursor, obj.cursor);
    }
}
