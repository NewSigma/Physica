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

namespace Physica::Core {
    template<ExprType type, class T1, class T2 = T1> class SparseVectorExpr;

    template<class VectorType, Scalar T>
    class SparseVectorExpr<ExprType::Mul, VectorType, T>
            : public RSparseVector<SparseVectorExpr<ExprType::Mul, VectorType, T>> {
        using This = SparseVectorExpr<ExprType::Mul, VectorType, T>;
        using Base = RSparseVector<This>;
        using typename Base::ScalarType;
        using typename Base::NonZeroPair;

        const VectorType& v;
        const T& s;
    public:
        SparseVectorExpr(const RSparseVector<VectorType>& v_, const T& s_)
                : v(v_.getDerived()), s(s_.getDerived()) {}
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
    template<class VectorType, Scalar T>
    [[nodiscard]] inline auto operator*(const RSparseVector<VectorType>& v, const T& x) noexcept {
        return SparseVectorExpr<ExprType::Mul, VectorType, T>(v.getDerived(), x);
    }

    template<class VectorType, Scalar T>
    [[nodiscard]] inline auto operator*(const T& x, const RSparseVector<VectorType>& v) noexcept {
        return SparseVectorExpr<ExprType::Mul, VectorType, T>(v * x);
    }
}

namespace Physica {
    using namespace Core;

    template<ExprType type, class Exp, Scalar T>
    class Traits<SparseVectorExpr<type, Exp, T>> {
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename Exp::ScalarType, T>::Type;
        constexpr static size_t SizeAtCompile = Exp::SizeAtCompile;
    };
}
