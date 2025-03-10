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
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"

namespace Physica {
    template<Matrix M, Vector U>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<M, U>> {
        using This = MatrixVectorProduct<M, U>;
        using Base = RValueVector<This>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const M& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const M& mat_, const U& vec_);
        MatrixVectorProduct(const This&) = default;
        MatrixVectorProduct(This&&) noexcept = default;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& target) const;

        [[nodiscard]] inline CoDiff<T> calc(size_t index) const;
        [[nodiscard]] inline Tv calc_value(size_t index) const;

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff);

        auto values() const noexcept { return mat.values() * vec.values(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const M& getLHS() const noexcept { return mat; }
        [[nodiscard]] const U& getRHS() const noexcept { return vec; }
    };

    template<Matrix M, Vector U>
    MatrixVectorProduct<M, U>::MatrixVectorProduct(const M& mat_, const U& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector U>
    template<Vector V, class Executor>
    inline void MatrixVectorProduct<M, U>::assign(V& target) const {
        if constexpr (isReverseDiff) {
            if constexpr (MatrixOption::isColMatrix<M>()) {
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
            if constexpr (MatrixOption::isColMatrix<M>()) {
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

    template<Matrix M, Vector U>
    inline auto MatrixVectorProduct<M, U>::calc(size_t index) const -> CoDiff<T> {
        return mat.row(index) * vec;
    }

    template<Matrix M, Vector U>
    inline auto MatrixVectorProduct<M, U>::calc_value(size_t index) const -> Tv {
        return mat.values().row(index) * vec.values();
    }

    template<Matrix M, Vector U>
    template<Vector V>
    void MatrixVectorProduct<M, U>::reverse(const V& grad_) const noexcept requires(isReverseDiff) {
        assert(grad_.getLength() == getLength());
        const auto& grad = grad_.values();
        if constexpr (ReverseDiff<M>) {
            if constexpr (MatrixOption::isRowMatrix<M>()) {
                for (size_t i = 0; i < mat.getRow(); ++i)
                    mat.row(i).reverse(grad[i] * vec.values());
            }
            else {
                for (size_t i = 0; i < mat.getCol(); ++i)
                    mat.col(i).reverse(grad * vec.calc_value(i));
            }
        }

        if constexpr (ReverseDiff<U>) {
            if constexpr (MatrixOption::isRowMatrix<M>()) {
                for (size_t i = 0; i < getLength(); ++i)
                    vec.reverse(mat.values().row(i) * grad.calc(i));
            }
            else
                vec.reverse(mat.values().transpose() * grad);
        }
    }

    template<Matrix M, Vector U>
    [[nodiscard]] inline auto operator*(const M& mat, const U& vec) noexcept requires(M::RowAtCompile != 1 && !CUDA<M> && !CUDA<U>) {
        return MatrixVectorProduct<M, U>(mat, vec);
    }

    template<Matrix M, Vector U>
    [[nodiscard]] inline auto operator*(const M& mat, const U& vec) requires(M::RowAtCompile == 1 && U::SizeAtCompile == 1 && !CUDA<M> && !CUDA<U>) {
        return mat.row(0) * vec;
    }
}

namespace Physica {
    template<Matrix M, Vector U>
    class Traits<MatrixVectorProduct<M, U>> {
        using M1 = std::remove_cvref_t<M>;
        using U1 = std::remove_cvref_t<U>;
        static_assert(M1::ColAtCompile == U1::SizeAtCompile || M1::ColAtCompile == Dynamic || U1::SizeAtCompile == Dynamic,
                "Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename U1::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile;
        constexpr static bool FastAssign = MatrixOption::isColMatrix<M>();
        constexpr static bool FastPacket = false;
    };
}
