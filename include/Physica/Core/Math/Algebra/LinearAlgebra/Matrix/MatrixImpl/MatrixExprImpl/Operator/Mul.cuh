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
    class device_obj<MatrixExpr<ExprID::Mul, M, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Mul, M, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Mul, M, U>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        using Base::operator*;
        [[nodiscard]] __host__ __device__ auto operator*(Scalar auto x) const noexcept;
        [[nodiscard]] __host__ __device__ auto operator-(this auto&&) noexcept;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto) const;

        [[nodiscard]] __host__ __device__ auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, U>>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (getRHS() * x);
    }

    template<Matrix M, Scalar U>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, U>>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * (-std::forward<Self>(self).getRHS());
    }

    template<Matrix M, Scalar U>
    __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, U>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col);
        else
            return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, U>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().value();
    }

    template<Matrix M, Vector V>
    class device_obj<MatrixExpr<ExprID::Mul, M, V>>
            : public device_obj<BinaryMatrixExpr<ExprID::Mul, M, V>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Mul, M, V>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] __host__ __device__ auto values(this auto&& self) noexcept;
    };

    template<Matrix M, Vector V>
    __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, V>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return Base::getLHS().calc(row, col, block) * Base::getRHS().calc(row, block);
    }

    template<Matrix M, Vector V>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Mul, M, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M1, Matrix M2>
    class device_obj<MatrixExpr<ExprID::Mul, M1, M2>>
            : public device_obj<BinaryMatrixExpr<ExprID::Mul, M1, M2>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Mul, M1, M2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] __host__ __device__ auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    __device__ auto device_obj<MatrixExpr<ExprID::Mul, M1, M2>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return getLHS().calc(row, col, block) * getRHS().calc(row, col, block);
    }

    template<Matrix M1, Matrix M2>
    void device_obj<MatrixExpr<ExprID::Mul, M1, M2>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(hadamard(rhs.values(), grad));
        if constexpr (Diffable<M2>)
            rhs.reverse(hadamard(lhs.values(), grad));
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto device_obj<MatrixExpr<ExprID::Mul, M1, M2>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator*(M&& m, U&& x) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Mul, remove_device_obj_t<M&&>, U&&>>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator*(U&& x, M&& m) noexcept requires(DeviceObj<M>) {
        return m * x;
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto hadamard(M&& m, V&& x) noexcept requires(DeviceObj<M> && DeviceObj<V>) {
        return device_obj<MatrixExpr<ExprID::Mul, remove_device_obj_t<M&&>, V&&>>(std::forward<M>(m), std::forward<V>(x));
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto hadamard(V&& x, M&& m) noexcept requires(DeviceObj<M> && DeviceObj<V>) {
        return hadamard(m, x);
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] auto hadamard(M1&& m1, M2&& m2) noexcept requires(DeviceObj<M1> && DeviceObj<M2>) {
        if constexpr (!canonicalized(m1, m2))
            return hadamard(std::forward<M2>(m2), std::forward<M1>(m1));
        else
            return device_obj<MatrixExpr<ExprID::Mul, remove_device_obj_t<M1&&>, remove_device_obj_t<M2&&>>>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
