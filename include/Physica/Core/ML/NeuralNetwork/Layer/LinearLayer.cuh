/*
 * Copyright 2024-2025 Weibo He.
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

#include "LinearLayer.h"
#include "LayerBase.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica {
    template<Scalar T, bool WithBias>
    class device_obj<LinearLayer<T, WithBias>> : public device_obj<LayerBase<LinearLayer<T, WithBias>>> {
        using host_obj = LinearLayer<T, WithBias>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<host_obj>>;

        using Tm = T::MachineType;
        using BiasType = std::conditional<WithBias, device_obj<VectorND<T>>, PlainStruct<void>>::type;
    public:
        template<Scalar U>
        using MatrixND = host_obj::template MatrixND<U>;
    private:
        using Tv = T::ValueType;

        device_obj<MatrixND<T>> weights;
        [[no_unique_address]] BiasType bias;
    public:
        device_obj() = default;
        device_obj(size_t inputDim, size_t outputDim);
        template<Scalar U>
        device_obj(const device_obj<LinearLayer<U, WithBias>>& layer);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Vector V>
        [[nodiscard]] inline CoDiff<device_obj<VectorND<T>>> forward(const V& x) const requires(CUDA<V>);
        template<Matrix M>
        [[nodiscard]] inline CoDiff<device_obj<MatrixND<T>>> forward(const M& x) const requires(CUDA<M>);
        void reverse(const This& __restrict other) const noexcept;

        void step(auto& optimizer);
        void zero_grad();

        void resize(size_t inputDim, size_t outputDim);
        [[nodiscard]] inline host_obj toHost() const;
        [[nodiscard]] inline host_obj toHostAsync() const;
        inline void toHost(host_obj& obj) const;
        inline void toHostAsync(host_obj& obj) const;

        void fill(Tv x);
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_xavier_uniform(Tv gain);
        template<RNG R>
        void random_xavier_normal(Tv gain);
        template<RNG R>
        void random_kaiming_uniform(Tv gain);
        template<RNG R>
        void random_kaiming_normal(Tv gain);
        template<RNG R>
        void random_any(auto& distribution);

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const auto& getWeights() const noexcept { return weights; }
        [[nodiscard]] auto& getWeights() noexcept { return weights; }
        [[nodiscard]] const auto& getBias() const noexcept requires(WithBias) { return bias; }
        [[nodiscard]] auto& getBias() noexcept requires(WithBias) { return bias; }
        /* Friends */
        friend class LinearLayer<T, WithBias>;
    };

    template<Scalar T, bool WithBias>
    device_obj<LinearLayer<T, WithBias>>::device_obj(size_t inputDim, size_t outputDim) {
        resize(inputDim, outputDim);
    }

    template<Scalar T, bool WithBias>
    template<Scalar U>
    device_obj<LinearLayer<T, WithBias>>::device_obj(const device_obj<LinearLayer<U, WithBias>>& layer) {
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
    template<Vector V>
    inline auto device_obj<LinearLayer<T, WithBias>>::forward(const V& x) const -> CoDiff<device_obj<VectorND<T>>> requires(CUDA<V>) {
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
    template<Matrix M>
    inline auto device_obj<LinearLayer<T, WithBias>>::forward(const M& x) const -> CoDiff<device_obj<MatrixND<T>>>  requires(CUDA<M>) {
        assert(x.getRow() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
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
    void device_obj<LinearLayer<T, WithBias>>::reverse(const This& __restrict other) const noexcept {
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
    void device_obj<LinearLayer<T, WithBias>>::step(auto& optimizer) {
        if constexpr (ReverseDiff<T>) {
            optimizer.step(weights);
            if constexpr (WithBias)
                optimizer.step(bias);
        }
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::zero_grad() {
        if constexpr (ReverseDiff<T>) {
            weights.grads().zeros();
            if constexpr (WithBias)
                bias.grads().zeros();
        }
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::resize(size_t inputDim, size_t outputDim) {
        weights.resize(outputDim, inputDim);
        if constexpr (WithBias)
            bias.resize(outputDim);
    }

    template<Scalar T, bool WithBias>
    inline auto device_obj<LinearLayer<T, WithBias>>::toHost() const -> host_obj {
        host_obj result = toHostAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, bool WithBias>
    inline auto device_obj<LinearLayer<T, WithBias>>::toHostAsync() const -> host_obj {
        host_obj result(getInputDim(), getOutputDim());
        toHostAsync(result);
        return result;
    }

    template<Scalar T, bool WithBias>
    inline void device_obj<LinearLayer<T, WithBias>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, bool WithBias>
    inline void device_obj<LinearLayer<T, WithBias>>::toHostAsync(host_obj& obj) const {
        weights.toHostAsync(obj.weights);
        if constexpr (WithBias)
            bias.toHostAsync(obj.bias);
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::fill(Tv x) {
        weights = x;
        if constexpr (WithBias)
            bias = x;
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_normal() {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_uniform(Tv gain) {
        const auto factor = (gain * sqrt(Tv(6) / Tv(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<Tm> dist(-factor, factor);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_normal(Tv gain) {
        const Tv deviation = gain * sqrt(Tv(2) / Tv(getInputDim() + getOutputDim()));
        random_normal<R>();
        weights *= deviation;
        if constexpr (WithBias)
            bias *= deviation;
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_kaiming_uniform(Tv gain) {
        const Tm bound = (gain * sqrt(Tv(3) / Tv(getInputDim()))).toMachine();
        std::uniform_real_distribution<Tm> dist(-bound, bound);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_kaiming_normal(Tv gain) {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();

        const Tv factor = gain / sqrt(Tv(getInputDim()));
        weights *= factor;
        if constexpr (WithBias)
            bias *= factor;
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_any(auto& distribution) {
        weights.random_any(distribution);
        if constexpr (WithBias)
            bias.random_any(distribution);
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, bool WithBias>
    const H5Group device_obj<LinearLayer<T, WithBias>>::read(const H5Loc& loc, const char* name) {
        host_obj obj{};
        auto group = obj.read(loc, name);
        obj.toDeviceAsync(*this);
        return group;
    }

    template<Scalar T, bool WithBias>
    H5Group device_obj<LinearLayer<T, WithBias>>::write(H5Loc& loc, const char* name) const {
        host_obj obj = toHost();
        return obj.write(loc, name);
    }
#endif

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }

    template<Scalar T, bool WithBias>
    auto LinearLayer<T, WithBias>::toDevice() const {
        auto result = toDeviceAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, bool WithBias>
    auto LinearLayer<T, WithBias>::toDeviceAsync() const {
        device_obj<This> result{};
        toDeviceAsync(result);
        return result;
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::toDevice(device_obj<This>& obj) const {
        toDeviceAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::toDeviceAsync(device_obj<This>& obj) const {
        weights.toDeviceAsync(obj.weights);
        if constexpr (WithBias)
            bias.toDeviceAsync(obj.bias);
    }
}

namespace Physica {
    template<Scalar T, bool WithBias>
    class Traits<device_obj<LinearLayer<T, WithBias>>> : public Traits<LinearLayer<T, WithBias>> {};
}
