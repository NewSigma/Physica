/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.cuh"

namespace Physica {
    template<Matrix M, Scalar U>
    class device_obj<MatrixExpr<ExprID::Add, M, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Add, M, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Add, M, U>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) + Base::getRHS();
        }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
    };

    template<Matrix M, Scalar U>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Add, M, U>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().value();
    }

    template<Matrix M, Vector U>
    class device_obj<MatrixExpr<ExprID::Add, M, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Add, M, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Add, M, U>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        __host__ __device__ void assign(Matrix auto& target) const {
            constexpr bool FastAssign = Traits<std::remove_cvref_t<M>>::FastAssign;
            if constexpr (FastAssign) {
                Base::getLHS().assign(target);
                target += Base::getRHS();
            }
            else
                Base::assign(target);
        }

        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) + Base::getRHS().calc(row);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row);
        }

        using Base::reverse;
        void reverse(const Matrix auto& grad) const noexcept {
            static_assert(isReverseDiff);
            if constexpr (ReverseDiff<M>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(grad.sum_cols());
        }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
    };

    template<Matrix M, Vector U>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Add, M, U>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M1, Matrix M2>
    class device_obj<MatrixExpr<ExprID::Add, M1, M2>>
            : public device_obj<BinaryMatrixExpr<ExprID::Add, M1, M2>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Add, M1, M2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) + Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row, col);
        }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Add, M1, M2>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(M&& m, U&& x) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Add, M&&, U&&>>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(U&& x, M&& m) noexcept requires(DeviceObj<M>) {
        return m + x;
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(M&& m, V&& x) noexcept requires(DeviceObj<M> && DeviceObj<V>) {
        return device_obj<MatrixExpr<ExprID::Add, M&&, V&&>>(std::forward<M>(m), std::forward<V>(x));
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(V&& x, M&& m) noexcept requires(DeviceObj<M> && DeviceObj<V>) {
        return m + x;
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(M1&& m1, M2&& m2) noexcept requires(DeviceObj<M1> && DeviceObj<M2>) {
        if constexpr (!canonicalized(m1, m2))
            return std::forward<M2>(m2) + std::forward<M1>(m1);
        else
            return device_obj<MatrixExpr<ExprID::Add, remove_device_obj_t<M1&&>, remove_device_obj_t<M2&&>>>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
