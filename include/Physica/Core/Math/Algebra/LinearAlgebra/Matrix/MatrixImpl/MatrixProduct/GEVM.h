/*
 * Copyright 2021-2026 Weibo He.
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
    template<Scalar T, size_t Order> class IdentityMatrix;
    /**
     * \class GEVM represents a outer product expression
     *
     * To compute (row vector) * (matrix), users should convert it to (matrix)^T * (col vector).
     */
    template<Vector V, Matrix M>
    class GEVM : public RValueMatrix<GEVM<V, M>> {
        using This = GEVM<V, M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> vec;
        LazyDestroy<M> mat;
    public:
        GEVM(V vec_, M mat_);
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

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
    };

    template<Vector V, Matrix M>
    GEVM<V, M>::GEVM(V vec_, M mat_) : vec(std::forward<V>(vec_)), mat(std::forward<M>(mat_)) {
        assert(vec.getLength() > 0 && mat.getCol() > 0 && "[Error]: Empty vector or matrix");
        assert(mat.getRow() == 1 && "[Error]: Dimensions do not match");
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign(Matrix auto& target) const {
        assign_base(target);
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_base(Matrix auto& target) const {
        if constexpr (MatrixMajor::isColMatrix<decltype(target)>()) {
            for (size_t c = 0; c < getCol(); ++c)
                target.col(c) = vec * mat(0, c);
        }
        else {
            for (size_t r = 0; r < getRow(); ++r)
                target.row(r) = vec.calc(r) * mat.flatten();
        }
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_add(Matrix auto& target) const {
        if constexpr (MatrixMajor::isColMatrix<decltype(target)>()) {
            for (size_t c = 0; c < getCol(); ++c)
                target.col(c) += vec * mat.calc(0, c);
        }
        else {
            for (size_t r = 0; r < getRow(); ++r)
                target.row(r) += vec.calc(r) * mat.flatten();
        }
    }

    template<Vector V, Matrix M>
    auto GEVM<V, M>::calc(size_t row, size_t col) const -> T {
        return vec.calc(row) * mat.calc(0, col);
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        if constexpr (Diffable<V>)
            vec.reverse(grad * mat.transpose());
        if constexpr (Diffable<M>)
            mat.reverse(grad.transpose() * vec);
    }

    template<Vector V, Matrix M>
    auto GEVM<V, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Vector V, Matrix M>
    auto&& GEVM<V, M>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Vector V, Matrix M>
    auto&& GEVM<V, M>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Vector V, Matrix M>
    __host__ __device__ consteval size_t GEVM<V, M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }

    template<Vector V, Matrix M>
    __host__ __device__ consteval size_t GEVM<V, M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }
}

namespace Physica {
    template<Vector V, Matrix M>
    class Traits<GEVM<V, M>> {
        using V1 = std::remove_cvref<V>::type;
        using M1 = std::remove_cvref<M>::type;
        using T1 = V1::ScalarType;
        using T2 = M1::ScalarType;

        static_assert(M1::getRowAtCompile() == 1 || M1::getRowAtCompile() == Dynamic, "[Error]: Outer product requires that the rows of M be 1");
        static_assert(!instanceof_tx<IdentityMatrix, M>, "[Error]: This pattern is unusual and is likely a bug. Consider rewriting if this is a false-positive");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        constexpr static int Major = MatrixMajor::BothMajor;
    };
}
