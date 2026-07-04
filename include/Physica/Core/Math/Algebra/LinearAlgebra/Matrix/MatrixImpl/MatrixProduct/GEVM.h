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
    class GEVM final : public RValueMatrix<GEVM<V, M>> {
        using This = GEVM<V, M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<V> vec;
        decay_rvalue_t<M> mat;
    public:
        GEVM(V vec_, M mat_);
        GEVM(const This&) = default;
        GEVM(This&&) noexcept = default;
        ~GEVM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        using Base::operator*;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_add(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] CoDiff<T> sum() const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
    private:
        void check_prod() const noexcept;
    };

    template<Vector V, Matrix M>
    GEVM<V, M>::GEVM(V vec_, M mat_) : vec(std::forward<V>(vec_)), mat(std::forward<M>(mat_)) {
        assert(!vec.empty() && mat.getCol() > 0 && "[Error]: Empty vector or matrix");
    }

    template<Vector V, Matrix M>
    auto GEVM<V, M>::operator*(this auto&& self, Vector auto&& v) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * (std::forward<Self>(self).getRHS() * std::forward<decltype(v)>(v));
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign(Matrix auto& target) const {
        assign_base(target);
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_base(Matrix auto& target) const {
        check_prod();
        if constexpr (MatrixMajor::isColMatrix<decltype(target)>()) {
            for (size_t c = 0; c < getCol(); ++c)
                target.col(c) = vec * mat.calc(0, c);
        }
        else {
            for (size_t r = 0; r < getRow(); ++r)
                target.row(r) = vec.calc(r) * mat.flatten();
        }
    }

    template<Vector V, Matrix M>
    void GEVM<V, M>::assign_add(Matrix auto& target) const {
        check_prod();
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
        check_prod();
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
    auto GEVM<V, M>::sum() const -> CoDiff<T> {
        return getLHS().sum() * getRHS().sum();
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

    template<Vector V, Matrix M>
    void GEVM<V, M>::check_prod() const noexcept {
        static_assert(!instanceof_tx<M, IdentityMatrix>, "[Error]: This pattern is unusual and is likely a bug. Consider drop it if this is a false-positive");
        if constexpr (mat.getRowAtCompile() == Dynamic)
            assert(mat.getRow() == 1 && "[Error]: Outer product requires that the rows of M be 1");
        else
            static_assert(mat.getRowAtCompile() == 1, "[Error]: Outer product requires that the rows of M be 1");
    }
}

namespace Physica {
    template<Vector V, Matrix M>
    class Traits<GEVM<V, M>> {
        using T1 = std::remove_cvref<V>::type::ScalarType;
        using T2 = std::remove_cvref<M>::type::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
    };
}
