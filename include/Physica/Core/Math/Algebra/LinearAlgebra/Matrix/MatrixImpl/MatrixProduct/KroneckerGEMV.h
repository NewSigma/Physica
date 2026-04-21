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

#include "Kronecker.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat, V&& vec);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto& target) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    GEMV<M, V>::GEMV(M&& mat, V&& vec) : mat(std::forward<M>(mat)), vec(std::forward<V>(vec)) {}

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const noexcept {
        target.assert_assign(*this);

        size_t row = mat.getRHS().getCol();
        size_t col = mat.getLHS().getCol();
        (mat.getRHS() * vec.reshape_col(row, col) * mat.getLHS().transpose()).flatten().assign(target);
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto& target) const noexcept {
        target.assert_assign(*this);

        size_t row = mat.getRHS().getCol();
        size_t col = mat.getLHS().getCol();
        (mat.getRHS() * vec.reshape_col(row, col) * mat.getLHS().transpose()).flatten().assign_add(target);
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    auto GEMV<M, V>::calc(size_t index) const -> CoDiff<T> {
        size_t row = mat.getRHS().getCol();
        size_t col = mat.getLHS().getCol();
        return (mat.getRHS() * vec.reshape_col(row, col) * mat.getLHS().transpose()).flatten().calc(index);
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    __host__ __device__ consteval bool GEMV<M, V>::isFastAssign() noexcept {
        return MatrixMajor::isColMatrix<M>();
    }

    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    __host__ __device__ consteval size_t GEMV<M, V>::getSizeAtCompile() noexcept {
        return std::max(std::remove_cvref_t<M>::RowAtCompile, std::remove_cvref_t<V>::getSizeAtCompile());
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<Kronecker, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::ColAtCompile == V1::getSizeAtCompile() || M1::ColAtCompile == Dynamic || V1::getSizeAtCompile() == Dynamic,
                "Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
    };
}
