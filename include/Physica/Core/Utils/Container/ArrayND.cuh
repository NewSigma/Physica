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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh" // IWYU pragma: export
#include "ArrayND.h"

namespace Physica {
    template<class T, int... Dims>
    class device_obj<ArrayND<T, Dims...>> {
        using host_obj = ArrayND<T, Dims...>;
        using This = device_obj<host_obj>;
        using ArrayType = device_obj<typename host_obj::ArrayType>;
        using IndexType = host_obj::IndexType;
        using ShapeType = host_obj::ShapeType;

        constexpr static bool StaticShape = host_obj::StaticShape;
    public:
        constexpr static int NDim = host_obj::NDim;
        constexpr static size_t SizeAtCompile = host_obj::SizeAtCompile;
        static_assert(NDim > 2, "[Error]: Invalid Dim");
    private:
        ArrayType arr;
        [[no_unique_address]] ShapeType shape;
    public:
        device_obj() = default;
        explicit __host__ __device__ device_obj(IndexType shape_, auto&&... args);
        explicit __host__ __device__ device_obj(std::integral auto... dims);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ auto& operator[](this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] __device__ auto& operator[](this auto&&, std::integral auto... dims) noexcept;
        /* Operations */
        __host__ __device__ void resize(IndexType shape_, auto&&... args);
        __host__ __device__ void resize(std::integral auto... dims);
        void reserve(size_t size);

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        [[nodiscard]] __host__ __device__ size_t toIndex1D(const IndexType& indices) const noexcept;
        [[nodiscard]] __host__ __device__ IndexType toIndexND(size_t index) const noexcept;

        __device__ void forND(std::invocable<T&, IndexType> auto func);
        __device__ void forND(std::invocable<const T&, IndexType> auto func) const;

        void zeros();
        void junk();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t dim(int index) const noexcept;
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data(this auto&&) noexcept;
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data_ptr(this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] __host__ __device__ auto&& asArray(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept;
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return arr.empty(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static int ndim() noexcept { return NDim; }
        /* Friends */
        friend class ArrayND<T, Dims...>;
    };
}

#include "ArrayImpl/ArrayNDImpl.cuh"
