/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<class MatrixType, class AnyScalar>
    class device_obj<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>>
            : public device_obj<RValueMatrix<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>>> {
        using host_obj = MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<MatrixType>& exp;
        const AnyScalar& scalar;
    public:
        __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) * ScalarType(scalar);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };

    template<class MatrixType, class ScalarType>
    [[nodiscard]] __device__ inline device_obj<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>>
    operator*(const ScalarBase<ScalarType>& s, const device_obj<RValueMatrix<MatrixType>>& m) noexcept {
        return {m, s};
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] __device__ inline device_obj<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>>
    operator*(const device_obj<RValueMatrix<MatrixType>>& m, const ScalarBase<ScalarType>& s) noexcept {
        return {m, s};
    }
}
