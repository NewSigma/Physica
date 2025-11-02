/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"
#include "RValueVector.h"

namespace Physica {
    template<class Derived>
    class device_obj<RValueVector<Derived>> : public CRTPBase<device_obj<RValueVector<Derived>>> {
        using host_obj = RValueVector<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<device_obj<Derived>>;
    public:
        using ScalarType = TraitsType::ScalarType;
        constexpr static size_t SizeAtCompile = TraitsType::SizeAtCompile;
        using PacketType = device_obj<BestPacket<ScalarType, SizeAtCompile>>::Type;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isDiffable = ScalarType::isDiffable;
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static int MaxThreadPerBlock = 256;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
    private:
        using RealsRtnTy = std::conditional<isComplex, device_obj<RealVector<Derived>>, device_obj<Derived>&>::type;
        using ValuesRtnTy = std::conditional<isDiffable, device_obj<ValueVector<Derived>>, device_obj<Derived>&>::type;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V>
        __host__ __device__ void assign(V& target) const requires(CUDA<V>);
        __host__ __device__ void assert_assign(const Vector auto& target) const noexcept;

        [[nodiscard]] __device__ T calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const { return Base::getDerived().calc_value(index); }
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept requires(isReverseDiff);

        __host__ __device__ void resize(const Vector auto& x) { resize(x.getLength()); }
        __host__ __device__ auto resize(size_t length) { return Base::getDerived().resize(length); }

        [[nodiscard]] __host__ __device__ auto transpose() const noexcept;

        [[nodiscard]] __device__ Tr norm() const;
        [[nodiscard]] __device__ Tr squaredNorm() const;
        [[nodiscard]] __device__ T max() const;
        [[nodiscard]] __device__ T min() const;
        [[nodiscard]] __host__ __device__ T sum() const;
        [[nodiscard]] __device__ T mean() const;
        [[nodiscard]] __device__ T variance() const;
        [[nodiscard]] __device__ T variance(const T& prior_mean) const;
        [[nodiscard]] __device__ T deviation() const;
        [[nodiscard]] __device__ T deviation(const T& prior_mean) const;
        [[nodiscard]] __device__ T lnSumExp() const;
        [[nodiscard]] __device__ T crossEntropy(size_t index) const;
        [[nodiscard]] __device__ T lnSoftmax(size_t index) const;
        [[nodiscard]] __device__ T softmax(size_t index) const;
        [[nodiscard]] __device__ auto crossProduct(const Vector auto& v) const noexcept requires(CUDA<decltype(v)>);

        [[nodiscard]] __device__ T max(int tid, int numThread, T* __restrict shared) const;
        [[nodiscard]] __device__ T sum(int tid, int numThread, T* __restrict shared) const;
        [[nodiscard]] __device__ T lnSumExp(int tid, int numThread, T* __restrict shared) const;
        [[nodiscard]] __device__ T mean(int tid, int numThread, T* __restrict shared) const;

        [[nodiscard]] __host__ __device__ RealsRtnTy reals() const noexcept;
        [[nodiscard]] __host__ __device__ auto imags() const noexcept;
        [[nodiscard]] __host__ __device__ auto squaredNorms() const noexcept;
        [[nodiscard]] __host__ __device__ auto norms() const noexcept;
        [[nodiscard]] __host__ __device__ ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grads() const noexcept;

        [[nodiscard]] __host__ __device__ KernelConfig makeKernelConfig() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ static KernelConfig makeKernelConfig(size_t length) noexcept;
        __host__ __device__ static void static_assert_assign(const Vector auto& target) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    private:
        template<Vector V>
        __device__ void assign_impl(V& target) const requires(CUDA<V>);
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RValueVector<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.cuh"
#include "RValueVectorImpl/CrossProduct.cuh"
#include "RValueVectorImpl/VectorConvert.cuh"
#include "RValueVectorImpl/InnerDot.cuh"
#include "VectorExpr.cuh"

