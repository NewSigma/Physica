/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    /**
     * \class GEVM represents a outer product expression
     *
     * To compute (row vector) * (matrix), users should convert it to (matrix)^T * (col vector).
     */
    template<Vector V, Matrix M>
    class GEVM : public RValueMatrix<GEVM<V, M>> {
        using This = GEVM<V, M>;
    public:
        using Base = RValueMatrix<This>;
        using Base::isReverseDiff;
        using typename Base::T;
    private:
        const V& vec;
        const M& mat;
    public:
        GEVM(const V& vec_, const M& mat_);
        GEVM(const This&) = default;
        GEVM(This&&) noexcept = default;
        ~GEVM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_add(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] const V& getLHS() const noexcept { return vec; }
        [[nodiscard]] const M& getRHS() const noexcept { return mat; }
    };

    template<Vector V, Matrix M>
    GEVM<V, M>::GEVM(const V& vec_, const M& mat_) : vec(vec_), mat(mat_) {
        assert(vec.getLength() > 0 && mat.getCol() > 0 && "[Error]: Empty vector or matrix");
        assert(mat.getRow() == 1 && "[Error]: Dimensions do not match");
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign(Matrix auto& target) const {
        assign_base(target);
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_base(Matrix auto& target) const {
        if constexpr (MatrixOption::isColMatrix<decltype(target)>()) {
            for (size_t c = 0; c < getCol(); ++c)
                target.col(c) = vec * mat(0, c);
        }
        else {
            for (size_t r = 0; r < getRow(); ++r)
                target.row(r) = vec[r] * mat.flatten();
        }
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_add(Matrix auto& target) const {
        if constexpr (MatrixOption::isColMatrix<decltype(target)>()) {
            for (size_t c = 0; c < getCol(); ++c)
                target.col(c) += vec * mat.calc(0, c);
        }
        else {
            for (size_t r = 0; r < getRow(); ++r)
                target.row(r) += vec[r] * mat.flatten();
        }
    }

    template<Vector V, Matrix M>
    auto GEVM<V, M>::calc(size_t row, size_t col) const -> T {
        return vec.calc(row) * mat.calc(0, col);
    }

    template<Vector V, Matrix M>
    [[nodiscard]] auto operator*(const V& vec, const M& mat) noexcept {
        static_assert(M::RowAtCompile == 1, "[Error]: Outer product requires that the rows of M be 1");
        return GEVM<V, M>(vec, mat);
    }
}

namespace Physica {
    template<Vector V, Matrix M>
    class Traits<GEVM<V, M>> {
        static_assert(M::RowAtCompile == 1, "[Error]: Outer product requires that the rows of M be 1");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename V::ScalarType, typename M::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = V::SizeAtCompile;
        constexpr static size_t ColAtCompile = M::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#include "MatrixProductImpl/GEVMExpr.h"
