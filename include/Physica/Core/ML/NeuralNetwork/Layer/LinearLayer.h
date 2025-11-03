/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h" // IWYU pragma: export
#include "Physica/PlainStruct.h"
#include "LayerBase.h"

namespace Physica {
    template<Scalar T, bool WithBias = true>
    class LinearLayer : public LayerBase<LinearLayer<T, WithBias>> {
        using This = LinearLayer<T, WithBias>;
        using Base = LayerBase<This>;

        using Tv = T::ValueType;
        using Tm = T::MachineType;
        using BiasType = std::conditional<WithBias, VectorND<T>, PlainStruct<void>>::type;
    public:
        template<Scalar U>
        using MatrixND = Base::template MatrixND<U>;
        using device_obj_type = device_obj<This>;
    private:
        MatrixND<T> weights;
        [[no_unique_address]] BiasType bias;
    public:
        LinearLayer() = default;
        LinearLayer(size_t inputDim, size_t outputDim);
        template<Scalar U>
        LinearLayer(const LinearLayer<U, WithBias>& layer);
        LinearLayer(const This&) = default;
        LinearLayer(This&&) noexcept = default;
        ~LinearLayer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] CoDiff<VectorND<T>> forward(const Vector auto& x) const;
        void reverse(const This& __restrict other) const noexcept;

        void step(auto& optimizer);
        void zero_grad();

        void resize(size_t inputDim, size_t outputDim);
        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        void fill(Tv x);
        template<RNG R = Random<>>
        void random_normal();
        template<RNG R = Random<>>
        void random_xavier_uniform(Tv gain);
        template<RNG R = Random<>>
        void random_xavier_normal(Tv gain);
        template<RNG R = Random<>>
        void random_kaiming_uniform(Tv gain);
        template<RNG R = Random<>>
        void random_kaiming_normal(Tv gain);
        template<RNG R = Random<>>
        void random_any(auto& distribution);

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(LinearLayer& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const auto& getWeights() const noexcept { return weights; }
        [[nodiscard]] const auto& getBias() const noexcept requires(WithBias) { return bias; }
    private:
        friend class device_obj<This>;
    };

    template<Scalar T, bool WithBias>
    LinearLayer<T, WithBias>::LinearLayer(size_t inputDim, size_t outputDim) {
        resize(inputDim, outputDim);
    }

    template<Scalar T, bool WithBias>
    template<Scalar U>
    LinearLayer<T, WithBias>::LinearLayer(const LinearLayer<U, WithBias>& layer) {
        if constexpr (Diffable<T>) {
            weights = layer.getWeights();
            if constexpr (WithBias)
                bias = layer.getBias();
        }
        else {
            weights = layer.getWeights().values();
            if constexpr (WithBias)
                bias = layer.getBias().values();
        }
    }

    template<Scalar T, bool WithBias>
    auto LinearLayer<T, WithBias>::forward(const Vector auto& x) const -> CoDiff<VectorND<T>> {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (ReverseDiff<T>) {
            auto expr1 = weights * x;
            if constexpr (WithBias) {
                auto expr2 = expr1 + bias;
                const auto result = co_yield expr2.values();
                expr2.reverse(result.grads());
            }
            else {
                const auto result = co_yield expr1.values();
                expr1.reverse(result.grads());
            }
        }
        else {
            if constexpr (WithBias)
                co_return weights * x + bias;
            else
                co_return weights * x;
        }
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::reverse(const This& __restrict other) const noexcept {
        assert(this != &other && "[Error]: Self reverse is invalid");
        assert(getInputDim() == other.getInputDim());
        assert(getOutputDim() == other.getOutputDim());
        if constexpr (ReverseDiff<T>) {
            weights.reverse(other.weights.grads());
            if constexpr (WithBias)
                bias.reverse(other.bias.grads());
        }
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::step(auto& optimizer) {
        if constexpr (ReverseDiff<T>) {
            optimizer.step(weights);
            if constexpr (WithBias)
                optimizer.step(bias);
        }
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::zero_grad() {
        if constexpr (ReverseDiff<T>) {
            weights.grads().zeros();
            if constexpr (WithBias)
                bias.grads().zeros();
        }
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::resize(size_t inputDim, size_t outputDim) {
        weights.resize(outputDim, inputDim);
        if constexpr (WithBias)
            bias.resize(outputDim);
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::fill(Tv x) {
        weights = x;
        if constexpr (WithBias)
            bias = x;
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_normal() {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_xavier_uniform(Tv gain) {
        const auto factor = (gain * sqrt(Tv(6) / Tv(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<Tm> dist(-factor, factor);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_xavier_normal(Tv gain) {
        const auto deviation = (gain * sqrt(Tv(2) / Tv(getInputDim() + getOutputDim()))).toMachine();
        std::normal_distribution<Tm> dist(0, deviation);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_kaiming_uniform(Tv gain) {
        const Tm bound = (gain * sqrt(Tv(3) / Tv(getInputDim()))).toMachine();
        std::uniform_real_distribution<Tm> dist(-bound, bound);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_kaiming_normal(Tv gain) {
        random_normal<R>();
        const Tv factor = gain / sqrt(Tv(getInputDim()));
        weights *= factor;
        if constexpr (WithBias)
            bias *= factor;
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void LinearLayer<T, WithBias>::random_any(auto& distribution) {
        weights.template random_any<R>(distribution);
        if constexpr (WithBias)
            bias.template random_any<R>(distribution);
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, bool WithBias>
    const H5Group LinearLayer<T, WithBias>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        weights.values().read(group, "Weight");
        if constexpr (Diffable<T>)
            weights.grads().resize(weights.values());

        if constexpr (WithBias) {
            bias.values().read(group, "Bias");
            if constexpr (Diffable<T>)
                bias.grads().resize(bias.values());
        }
        return group;
    }

    template<Scalar T, bool WithBias>
    H5Group LinearLayer<T, WithBias>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        weights.values().write(group, "Weight");
        if constexpr (WithBias)
            bias.values().write(group, "Bias");
        return group;
    }
#endif

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::swap(LinearLayer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }
}

namespace Physica {
    template<Scalar T, bool B>
    class Traits<LinearLayer<T, B>> {
    public:
        using ScalarType = T;
        constexpr static bool WithBias = B;
    };
}
