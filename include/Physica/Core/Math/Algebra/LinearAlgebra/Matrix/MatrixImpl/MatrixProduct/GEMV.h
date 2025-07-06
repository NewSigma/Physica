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
#include "Physica/Core/Parallel/Parallel.h"

namespace Physica {
    template<Matrix M, Vector U>
    class GEMV : public RValueVector<GEMV<M, U>> {
        using This = GEMV<M, U>;
        using Base = RValueVector<This>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const LazyDestroy<M> mat;
        const LazyDestroy<U> vec;
    public:
        GEMV(M&& mat_, U&& vec_);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        inline void assign(Vector auto& target) const;

        [[nodiscard]] inline CoDiff<T> calc(size_t index) const;
        [[nodiscard]] inline Tv calc_value(size_t index) const;

        using Base::reverse;
        void reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff);

        auto values() const noexcept { return mat.values() * vec.values(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Matrix M, Vector U>
    GEMV<M, U>::GEMV(M&& mat_, U&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<U>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector U>
    template<ExecutePolicy P>
    inline void GEMV<M, U>::assign(Vector auto& target) const {
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
    inline auto GEMV<M, U>::calc(size_t index) const -> CoDiff<T> {
        return mat.row(index) * vec;
    }

    template<Matrix M, Vector U>
    inline auto GEMV<M, U>::calc_value(size_t index) const -> Tv {
        return mat.values().row(index) * vec.values();
    }

    template<Matrix M, Vector U>
    void GEMV<M, U>::reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff) {
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

    template<class Derived>
    template<Vector V>
    auto RValueMatrix<Derived>::operator*(V&& v) const& noexcept requires(RowAtCompile != 1 && !CUDA<V>) {
        return GEMV<const Derived&, V&&>(Base::getDerived(), std::forward<V>(v));
    }

    template<class Derived>
    template<Vector V>
    auto RValueMatrix<Derived>::operator*(V&& v) && noexcept requires(RowAtCompile != 1 && !CUDA<V>) {
        return GEMV<Derived&&, V&&>(std::move(Base::getDerived()), std::forward<V>(v));
    }

    template<class Derived>
    template<Vector V>
    auto RValueMatrix<Derived>::operator*(const V& v) const noexcept requires(RowAtCompile == 1 && !CUDA<V>) {
        return row(0) * v;
    }
}

namespace Physica {
    template<Matrix M, Vector U>
    class Traits<GEMV<M, U>> {
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
