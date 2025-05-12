/*
 * Copyright 2025 WeiBo He.
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

#include "DQMC.h"

namespace Physica {
    template<Scalar T>
    class GreenSampler {
        using This = GreenSampler<T>;
        using Tf = Diff<T, DiffMode::Forward, 1>;
    public:
        enum Observable {
            Density,
            DoubleOccupy
        };
    private:
        VectorND<Tf> samples;
        VectorND<T> signs;
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
        void sample(const DQMC<T>& dqmc, Observable type);
        void sample(const DQMC<T>& dqmc, T observe);
        [[nodiscard]] T calc();
        [[nodiscard]] T calcSign() const;

        void reset() { cursor = 0; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getSamples() const noexcept { return samples; }
        [[nodiscard]] const auto& getSigns() const noexcept { return signs; }
        [[nodiscard]] size_t getNumSample() const noexcept { return samples.getLength(); }
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(size_t numSample) : samples(numSample), signs(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    void GreenSampler<T>::sample(const DQMC<T>& dqmc, Observable type) {
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
        sample(dqmc, observe);
    }

    template<Scalar T>
    void GreenSampler<T>::sample(const DQMC<T>& dqmc, T observe) {
        samples[cursor].value() = dqmc.getLnPartitionZ();
        samples[cursor].grad() = observe;
        signs[cursor] = dqmc.getSign();
        cursor = (cursor + 1) % getNumSample();
    }

    template<Scalar T>
    T GreenSampler<T>::calc() {
        samples -= samples.values().max();
        const Tf sum = exp(samples) * signs;
        return sum.grad() / sum.value(); // Avoid ln, its grad is well defined but value may give NAN
    }

    template<Scalar T>
    T GreenSampler<T>::calcSign() const {
        VectorND<Tf> copy = samples;
        copy.grads() = signs;
        return copy.lnSumExp().grad();
    }

    template<Scalar T>
    void GreenSampler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        signs.swap(obj.signs);
        std::swap(cursor, obj.cursor);
    }
}
