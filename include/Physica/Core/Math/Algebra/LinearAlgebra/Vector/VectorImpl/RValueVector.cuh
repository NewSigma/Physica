/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Parallel/ThreadBlock.cuh"
#include "Physica/Core/Utils/CUDA/device_obj.h"
#include "RValueVector.h"

namespace Physica {
    template<class Derived, Scalar ScalarT>
    class device_obj<RValueVector<Derived, ScalarT>> : public CRTPBase<device_obj<RValueVector<Derived, ScalarT>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is not allowed");
        using host_obj = RValueVector<Derived, ScalarT>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = ScalarT;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator*(this auto&&, Scalar auto&& x) noexcept;
        [[nodiscard]] __host__ __device__ auto operator*(this auto&&, Matrix auto&& m) noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(this auto&&) noexcept;
        /* Operations */
        __host__ __device__ void assign(this const auto& self, Vector auto&& target);
        __device__ void assign(this const auto&, Vector auto&& target, instanceof_x<ThreadBlock> auto block);
        __host__ __device__ void assign_add(Vector auto& target) const;
        __host__ __device__ void assign_add_base(Vector auto& target) const;
        __host__ __device__ void assert_assign(const Vector auto& source) const noexcept;
        [[nodiscard]] __host__ __device__ constexpr KernelConfig makeKernelConfig() const noexcept;

        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index, instanceof_x<ThreadBlock> auto block) const;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(this const auto&, size_t index) noexcept;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(this const auto&, size_t index, size_t count) noexcept;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;

        __host__ __device__ void resize(const Vector auto& x) { resize(x.getLength()); }
        __host__ __device__ auto resize(size_t length) { return Base::getDerived().resize(length); }

        [[nodiscard]] __host__ __device__ decltype(auto) conjugate(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto transpose(this auto&&) noexcept;

        [[nodiscard]] __device__ Tr norm() const;
        [[nodiscard]] __device__ Tr squaredNorm(this const auto&);
        [[nodiscard]] __device__ T max() const;
        [[nodiscard]] __device__ T min() const;
        [[nodiscard]] __host__ __device__ T sum() const;
        [[nodiscard]] __device__ T mean() const;
        [[nodiscard]] __device__ T variance() const;
        [[nodiscard]] __device__ T variance(const T& prior_mean) const;
        [[nodiscard]] __device__ T deviation() const;
        [[nodiscard]] __device__ T deviation(const T& prior_mean) const;
        [[nodiscard]] __device__ T lnSumExp() const;
        [[nodiscard]] __device__ T crossEntropy(this const auto&, size_t index);
        [[nodiscard]] __device__ T lnSoftmax(size_t index) const;
        [[nodiscard]] __device__ T softmax(size_t index) const;
        [[nodiscard]] __host__ __device__ T prod(this const auto&) noexcept;

        [[nodiscard]] __device__ T max(this const auto&, instanceof_x<ThreadBlock> auto block);
        [[nodiscard]] __device__ T min(this const auto&, instanceof_x<ThreadBlock> auto block);
        [[nodiscard]] __device__ T sum(this const auto&, instanceof_x<ThreadBlock> auto block);
        [[nodiscard]] __device__ T mean(instanceof_x<ThreadBlock> auto block) const;
        [[nodiscard]] __device__ T lnSumExp(instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] __host__ __device__ decltype(auto) reals(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto imags(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto squaredNorms(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto norms(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ decltype(auto) grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isLValueVector() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isSparse() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastPacket() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile(const Vector auto& hint) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Scalar auto& source) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Vector auto& source) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T, Scalar S>
    class Traits<device_obj<RValueVector<T, S>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.cuh"
#include "RValueVectorImpl/Conjugate.cuh"
#include "RValueVectorImpl/Cross.cuh"
#include "RValueVectorImpl/Dot.cuh"
#include "RValueVectorImpl/VectorConvert/RealVector.cuh"
#include "RValueVectorImpl/VectorConvert/VectorConvert.cuh"
#include "VectorExpr.cuh"

