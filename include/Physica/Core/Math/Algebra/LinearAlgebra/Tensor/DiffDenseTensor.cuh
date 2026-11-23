/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Utils/Container/ArrayND.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/CompactTensor.cuh"
#include "DenseTensor.cuh"
#include "DiffDenseTensor.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    class device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>
            : public device_obj<CompactTensor<DenseTensor<Diff<T, Mode, Order>, Dims...>>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>>, Empty>::type {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
        using host_obj = DenseTensor<Diff<T, Mode, Order>, Dims...>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactTensor<host_obj>>;
    public:
        using typename Base::IndexType;
        using typename Base::ScalarType;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
        template<Scalar U>
        using rebind_scalar = device_obj<DenseTensor<U, Dims...>>;
    private:
        using ValueTensor = device_obj<DenseTensor<T, Dims...>>;
        using GradType = ScalarType::GradType;
        using GradTensorType = std::conditional<Order == 1, ValueTensor, device_obj<DenseTensor<GradType, Dims...>>>::type;

        ValueTensor v;
        GradTensorType g;
    public:
        device_obj() = default;
        __host__ __device__ explicit device_obj(IndexType shape);
        __host__ __device__ device_obj(IndexType shape, T init);
        __host__ __device__ device_obj(ValueTensor v_, GradTensorType g_);
        template<Tensor X>
        explicit(isReverseDiff()) device_obj(const X& x) requires(!ReverseDiff<X>);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(this auto&, IndexType shape);
        void reserve(this auto&, size_t size) noexcept;

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        template<RNG R>
        void random_uniform(this auto&);
        template<RNG R>
        void random_normal(this auto&);
        template<RNG R>
        void random_any(this auto&, auto& distribution);

        void zero_grad(this auto&);
        void junk(this auto&);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ constexpr IndexType getShape() const noexcept { return v.getShape(); }
        [[nodiscard]] __host__ __device__ constexpr size_t dim(int index) const noexcept { return v.dim(index); }
        [[nodiscard]] __host__ __device__ constexpr size_t getSize() const noexcept { return v.getSize(); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return v.empty(); }

        [[nodiscard]] __host__ __device__ auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto&& grads(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticShape() noexcept { return std::remove_cvref_t<ValueTensor>::isStaticShape(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<ValueTensor>::getSizeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static auto getShapeAtCompile() noexcept -> IndexType { return std::remove_cvref_t<ValueTensor>::getShapeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static auto getStrideAtCompile() noexcept -> IndexType { return std::remove_cvref_t<ValueTensor>::getStrideAtCompile(); }

        template<RNG R>
        [[nodiscard]] static This random_uniform(IndexType shape);
        template<RNG R>
        [[nodiscard]] static This random_normal(IndexType shape);
        template<RNG R>
        [[nodiscard]] static This random_any(IndexType shape, auto& distribution);
        /* Friends */
        friend class DenseTensor<Diff<T, Mode, Order>, Dims...>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    class Traits<device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>>
            : public Traits<DenseTensor<Diff<T, Mode, Order>, Dims...>> {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
    };
}

#include "TensorImpl/DiffDenseTensorImpl.cuh"
