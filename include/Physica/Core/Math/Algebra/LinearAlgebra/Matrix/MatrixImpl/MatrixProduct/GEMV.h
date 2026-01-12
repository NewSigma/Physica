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

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M, Vector V>
    class GEMV : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M mat_, V vec_);
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
        [[nodiscard]] Tv calc_value(size_t index) const;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept { return mat.values() * vec.values(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Matrix M, Vector V>
    GEMV<M, V>::GEMV(M mat_, V vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {}

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const noexcept {
        if constexpr (MatrixOption::isColMatrix<M>()) {
            size_t length = vec.getLength();
            (mat.col(0) * vec.calc(0)).template assign<P>(target);
            for (size_t i = 1; i < length; ++i)
                (mat.col(i) * vec.calc(i)).template assign_add<P>(target);
        }
        else {
            size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                target[i] = calc(i);
        }
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto& target) const noexcept {
        if constexpr (MatrixOption::isColMatrix<M>()) {
            size_t length = vec.getLength();
            for (size_t i = 0; i < length; ++i)
                (mat.col(i) * vec.calc(i)).template assign_add<P>(target);
        }
        else {
            size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                target[i] += calc(i);
        }
    }

    template<Matrix M, Vector V>
    auto GEMV<M, V>::calc(size_t index) const -> CoDiff<T> {
        return mat.row(index) * vec;
    }

    template<Matrix M, Vector V>
    auto GEMV<M, V>::calc_value(size_t index) const -> Tv {
        return mat.values().row(index) * vec.values();
    }

    template<Matrix M, Vector V>
    void GEMV<M, V>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        vec.assert_assign(grad);
        const auto& g = grad.values();
        if constexpr (ReverseDiff<M>) {
            if constexpr (MatrixOption::isRowMatrix<M>()) {
                for (size_t i = 0; i < mat.getRow(); ++i)
                    mat.row(i).reverse(g[i] * vec.values());
            }
            else {
                for (size_t i = 0; i < mat.getCol(); ++i)
                    mat.col(i).reverse(g * vec.calc_value(i));
            }
        }

        if constexpr (ReverseDiff<V>) {
            if constexpr (MatrixOption::isRowMatrix<M>()) {
                for (size_t i = 0; i < getLength(); ++i)
                    vec.reverse(mat.values().row(i) * g.calc(i));
            }
            else
                vec.reverse(mat.values().transpose() * g);
        }
    }

    template<class Derived>
    template<Vector V>
    [[nodiscard]] auto RValueMatrix<Derived>::operator*(this auto&& self, V&& v) noexcept requires(RowAtCompile != 1 && !CUDA<V>) {
        assert(self.getCol() == v.getLength());
        using Self = decltype(self);
        return GEMV<Self, V&&>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<class Derived>
    template<Vector V>
    [[nodiscard]] auto RValueMatrix<Derived>::operator*(const V& v) const noexcept requires(RowAtCompile == 1 && !CUDA<V>) {
        assert(getCol() == v.getLength());
        return row(0) * v;
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::ColAtCompile == V1::SizeAtCompile || M1::ColAtCompile == Dynamic || V1::SizeAtCompile == Dynamic,
                "Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile;
        constexpr static bool FastAssign = MatrixOption::isColMatrix<M>();
        constexpr static bool FastPacket = false;
    };
}

#include "GEMVExpr.h"
