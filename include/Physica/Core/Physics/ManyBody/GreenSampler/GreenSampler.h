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
    private:
        VectorND<T> lnAbsDets;
        VectorND<Tv> signs;
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
        [[nodiscard]] T calcMean(const VectorND<T>& observes, const VectorND<T>& lnWeights) const;
        [[nodiscard]] Tv calcSign() const noexcept;
        [[nodiscard]] Tv calcSign(const VectorND<T>& lnWeights) const;
        [[nodiscard]] T lnPartitionZ() const;

        void reset() { cursor = 0; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getLnAbsDets() const noexcept { return lnAbsDets; }
        [[nodiscard]] const auto& getSigns() const noexcept { return signs; }
        [[nodiscard]] size_t getNumSample() const noexcept { return signs.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    protected:
        void sample(T lnAbsDet, T sign) noexcept;
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(size_t numSample) : lnAbsDets(numSample), signs(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    T GreenSampler<T>::calcMean(const VectorND<T>& observes) const {
        return hadamard(observes, signs).mean() / calcSign();
    }

    template<Scalar T>
    T GreenSampler<T>::calcMean(const VectorND<T>& observes, const VectorND<T>& lnWeights) const {
        VectorND<T> buffer = hadamard(exp(lnWeights - lnWeights.max()), signs);
        return hadamard(buffer, observes).mean() / buffer.mean();
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign() const noexcept -> Tv {
        return signs.mean();
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign(const VectorND<T>& lnWeights) const -> Tv {
        VectorND<T> buffer = exp(lnWeights - lnWeights.max());
        return hadamard(buffer, signs).mean() / buffer.mean();
    }

    template<Scalar T>
    T GreenSampler<T>::lnPartitionZ() const {
        const T factor = lnAbsDets.max();
        return ln(exp(lnAbsDets - factor) * getSigns()) + factor;
    }

    template<Scalar T>
    void GreenSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lnAbsDets.swap(obj.lnAbsDets);
        signs.swap(obj.signs);
        std::swap(cursor, obj.cursor);
    }

    template<Scalar T>
    void GreenSampler<T>::sample(T lnAbsDet, T sign) noexcept {
        lnAbsDets[cursor] = lnAbsDet;
        signs[cursor] = sign;
        cursor = (cursor + 1) % getNumSample();
    }
}
