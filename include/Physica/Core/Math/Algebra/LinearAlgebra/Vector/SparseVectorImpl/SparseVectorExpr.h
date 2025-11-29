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

#include "RSparseVector.h"

namespace Physica {
    template<ExprID, Vector T1, class T2 = T1> class SparseVectorExpr;

    template<Vector T, Scalar U>
    class SparseVectorExpr<ExprID::Mul, T, U>
            : public RSparseVector<SparseVectorExpr<ExprID::Mul, T, U>> {
        using This = SparseVectorExpr<ExprID::Mul, T, U>;
        using Base = RSparseVector<This>;
        using typename Base::ScalarType;
        using typename Base::NonZeroPair;

        const T& v;
        const U& s;
    public:
        SparseVectorExpr(const RSparseVector<T>& v_, const U& s_) : v(v_.getDerived()), s(s_) {}
        SparseVectorExpr(const This&) = delete;
        SparseVectorExpr(This&&) noexcept = delete;
        ~SparseVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return v.calc(index) * s; }
        [[nodiscard]] NonZeroPair calcNonZero(size_t index) const {
            const auto pair = v.calcNonZero(index);
            return std::make_pair(pair.first, ScalarType(pair.second * s));
        }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] size_t getNumNonZero() const noexcept { return v.getNumNonZero(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Vector T, Scalar U>
    [[nodiscard]] auto operator*(const RSparseVector<T>& v, const U& x) noexcept {
        return SparseVectorExpr<ExprID::Mul, T, U>(v.getDerived(), x);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] auto operator*(const U& x, const RSparseVector<T>& v) noexcept {
        return SparseVectorExpr<ExprID::Mul, T, U>(v * x);
    }
}

namespace Physica {
    template<ExprID ID, Vector T, Scalar U>
    class Traits<SparseVectorExpr<ID, T, U>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, U>::Type;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
