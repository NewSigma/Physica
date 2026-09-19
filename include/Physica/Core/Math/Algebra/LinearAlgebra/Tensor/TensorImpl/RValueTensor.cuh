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

#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Parallel/ThreadBlock.cuh"
#include "Physica/Core/Utils/CUDA/device_obj.h"
#include "RValueTensor.h"

namespace Physica {
    template<class Derived, Scalar ScalarT>
    class device_obj<RValueTensor<Derived, ScalarT>> : public CRTP<device_obj<RValueTensor<Derived, ScalarT>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is not allowed");
        using host_obj = RValueTensor<Derived, ScalarT>;
        using This = device_obj<host_obj>;
        using Base = CRTP<This>;
    public:
        constexpr static int NDim = host_obj::NDim;
        using ScalarType = ScalarT;
        using IndexType = host_obj::IndexType;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void assign(this const auto& self, Tensor auto&& target);
        __host__ __device__ void assert_assign(const Tensor auto& source) const noexcept;
        [[nodiscard]] __host__ __device__ constexpr KernelConfig makeKernelConfig() const noexcept;

        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& indices) const;
        [[nodiscard]] __device__ decltype(auto) calc(std::integral auto... dims) const;
        [[nodiscard]] __host__ __device__ auto toIndex1D(const IndexType& indices) const noexcept;
        [[nodiscard]] __host__ __device__ auto toIndexND(size_t index) const noexcept;

        __host__ __device__ void resize(this auto&, const Tensor auto& x);
        __host__ __device__ auto resize(this auto&, std::integral auto... dims);
        __host__ __device__ auto resize(this auto&, IndexType shape);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept;
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int ndim() noexcept { return NDim; }
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isLValueTensor() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isSparse() noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T, Scalar S>
    class Traits<device_obj<RValueTensor<T, S>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueTensorImpl/RValueTensorImpl.cuh"
#include "TensorExpr.cuh" // IWYU pragma: export
