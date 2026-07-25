/*
 * Copyright 2023-2026 Weibo He.
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

#include "../CompactMatrix.cuh"
#include "Transpose.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class device_obj<Transpose<M>> : public device_obj<CompactMatrix<Transpose<M>>> {
        using host_obj = Transpose<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return transpose().getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return transpose().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return transpose().getOrder(); }
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto&& device_obj<Transpose<M>>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto device_obj<Transpose<M>>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values().transpose();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto device_obj<Transpose<M>>::data(this auto&& self) noexcept {
        return self.mat.getDerived().data();
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    class device_obj<Transpose<V>> : public device_obj<CompactMatrix<Transpose<V>>> {
        using host_obj = Transpose<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactMatrix<host_obj>>;
        using Ref = add_device_obj<V>::type;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
    public:
        __host__ __device__ explicit device_obj(Ref vec_) : vec(asStruct(vec_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    };

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ auto&& device_obj<Transpose<V>>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.vec.getDerived());
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ auto device_obj<Transpose<V>>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values();
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ auto device_obj<Transpose<V>>::data(this auto&& self) noexcept {
        return self.vec.getDerived().data();
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ size_t device_obj<Transpose<V>>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return 1;
    }
}
