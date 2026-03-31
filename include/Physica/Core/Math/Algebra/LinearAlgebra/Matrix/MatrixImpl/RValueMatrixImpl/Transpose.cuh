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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<class T>
    struct remove_transpose<device_obj<Transpose<T>>> {
        using Type = device_obj<T>;
    };

    template<Matrix M>
    class device_obj<Transpose<M>> : public device_obj<RValueMatrix<Transpose<M>>> {
        using host_obj = Transpose<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;

        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        __host__ __device__ device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return getExpr().calc(col, row); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return getExpr().calc_value(col, row); }

        void reverse(const Vector auto& grad) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getExpr().getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getExpr().getRow(); }
    };

    template<Matrix M>
    void device_obj<Transpose<M>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        getExpr().reverse(grad.transpose());
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<Transpose<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat.getDerived());
    }

    template<Vector V>
    class device_obj<Transpose<V>> : public device_obj<RValueMatrix<Transpose<V>>> {
        using host_obj = Transpose<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<V>::type;

        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        __host__ __device__ explicit device_obj(Ref vec_) : vec(asStruct(vec_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.getDerived().calc(col); }
        [[nodiscard]] __device__ Tv calc_value([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc_value(col); }

        void reverse(const Vector auto& grad) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getDerived().getLength(); }
    };

    template<Vector V>
    void device_obj<Transpose<V>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        vec.getDerived().reverse(grad);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<Transpose<M>>> : public Traits<Transpose<M>> {};

    template<Vector V>
    class Traits<device_obj<Transpose<V>>> : public Traits<Transpose<V>> {};
}
