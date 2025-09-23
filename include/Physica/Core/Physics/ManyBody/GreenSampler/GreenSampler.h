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

        static_assert(!T::isComplex, "[Error]: Observables are real");
    protected:
        using Tv = T::ValueType;
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
        [[nodiscard]] T calcMean(const VectorND<T>& observes) const;
        [[nodiscard]] T calcMeanWeighted(const VectorND<T>& observes) const;
        [[nodiscard]] Tv calcSign() const noexcept;
        [[nodiscard]] Tv calcSignWeighted() const;
        [[nodiscard]] T lnPartitionZ() const;

        void reset() { cursor = 0; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getLnPartitionZs() const noexcept { return samples.values(); }
        [[nodiscard]] const auto& getSigns() const noexcept { return samples.grads(); }
        [[nodiscard]] size_t getNumSample() const noexcept { return samples.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    protected:
        void sample(T lnZ, T sign) noexcept;
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(size_t numSample) : samples(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    T GreenSampler<T>::calcMean(const VectorND<T>& observes) const {
        return hadamard(observes, getSigns()).mean() / calcSign();
    }

    template<Scalar T>
    T GreenSampler<T>::calcMeanWeighted(const VectorND<T>& observes) const {
        const T factor = getLnPartitionZs().max();
        VectorND<Tf> buffer(getLnPartitionZs() - factor, observes);
        const Tf sum = exp(buffer) * getSigns();
        return sum.grad() / sum.value(); // Avoid ln, its grad is well defined but value may give NAN
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign() const noexcept -> Tv {
        return getSigns().mean();
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSignWeighted() const -> Tv {
        return samples.lnSumExp().grad();
    }

    template<Scalar T>
    T GreenSampler<T>::lnPartitionZ() const {
        const T factor = getLnPartitionZs().max();
        return ln(exp(getLnPartitionZs() - factor) * getSigns()) + factor;
    }

    template<Scalar T>
    void GreenSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        std::swap(cursor, obj.cursor);
    }

    template<Scalar T>
    void GreenSampler<T>::sample(T lnZ, T sign) noexcept {
        samples.values()[cursor] = lnZ;
        samples.grads()[cursor] = sign;
        cursor = (cursor + 1) % getNumSample();
    }
}
