/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    template<Matrix T, Vector U>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<T, U>> {
        using This = MatrixVectorProduct<T, U>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    private:
        const T& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const T& mat_, const U& vec_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept requires(ReverseDiff<ScalarType>) = default;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V, class Executor = SequentialExecutor>
        inline void assignTo(V& target) const;

        [[nodiscard]] inline CoDiff<ScalarType> calc(size_t index) const;
        [[nodiscard]] inline ValueType calc_value(size_t index) const;

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const T& getLHS() const noexcept { return mat; }
        [[nodiscard]] const U& getRHS() const noexcept { return vec; }
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<T, U>::MatrixVectorProduct(const T& mat_, const U& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix T, Vector U>
    template<Vector V, class Executor>
    inline void MatrixVectorProduct<T, U>::assignTo(V& target) const {
        if constexpr (isReverseDiff) {
            if constexpr (MatrixOption::isColMatrix<T>()) {
                target = mat.values().col(0) * vec.values().calc(0);
                for (size_t i = 1; i < vec.getLength(); ++i)
                    target += mat.values().col(i) * vec.values().calc(i);
            }
            else {
                for (size_t i = 0; i < getLength(); ++i)
                    target[i] = calc_value(i);
            }
        }
        else {
            if constexpr (MatrixOption::isColMatrix<T>()) {
                target = mat.col(0) * vec.calc(0);
                for (size_t i = 1; i < vec.getLength(); ++i)
                    target += mat.col(i) * vec.calc(i);
            }
            else {
                for (size_t i = 0; i < getLength(); ++i)
                    target[i] = calc(i);
            }
        }
    }

    template<Matrix T, Vector U>
    inline auto MatrixVectorProduct<T, U>::calc(size_t index) const -> CoDiff<ScalarType> {
        return mat.row(index) * vec;
    }

    template<Matrix T, Vector U>
    inline auto MatrixVectorProduct<T, U>::calc_value(size_t index) const -> ValueType {
        return mat.values().row(index) * vec.values();
    }

    template<Matrix T, Vector U>
    template<Vector V>
    void MatrixVectorProduct<T, U>::reverse(const V& grad_) const noexcept requires(isReverseDiff) {
        const auto& grad = grad_.values();
        if constexpr (ReverseDiff<T>)
            mat.reverse(grad * vec.values().transpose());
        if constexpr (ReverseDiff<U>)
            vec.reverse(mat.values().transpose() * grad);
    }

    template<Matrix T, Vector U>
    [[nodiscard]] inline auto operator*(const T& mat, const U& vec) noexcept requires(T::RowAtCompile != 1) {
        static_assert(T::ColAtCompile == U::SizeAtCompile ||
                      T::ColAtCompile == Dynamic ||
                      U::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
        return MatrixVectorProduct<T, U>(mat, vec);
    }

    template<Matrix T, Vector U>
    [[nodiscard]] inline auto operator*(const T& mat, const U& vec) requires(T::RowAtCompile == 1 && U::SizeAtCompile == 1) {
        return mat.row(0) * vec;
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<MatrixVectorProduct<T, U>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = T::RowAtCompile;
        constexpr static bool FastAssign = MatrixOption::isColMatrix<T>();
        constexpr static bool FastPacket = false;
    };
}
