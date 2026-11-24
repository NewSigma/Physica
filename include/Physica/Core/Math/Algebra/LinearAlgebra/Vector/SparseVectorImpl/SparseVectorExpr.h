/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "RSparseVector.h"

namespace Physica {
    template<ExprID, Vector T1, class T2 = T1> class SparseVectorExpr;

    template<Vector V, Scalar U>
    class SparseVectorExpr<ExprID::Mul, V, U>
            : public RSparseVector<SparseVectorExpr<ExprID::Mul, V, U>> {
        using This = SparseVectorExpr<ExprID::Mul, V, U>;
        using Base = RSparseVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::NonZeroPair;
    private:
        const V& v;
        const U& s;
    public:
        SparseVectorExpr(const RSparseVector<V>& v_, const U& s_) : v(v_.getDerived()), s(s_) {}
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
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] constexpr size_t getNumNonZero() const noexcept { return v.getNumNonZero(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Vector V, Scalar U>
    consteval size_t SparseVectorExpr<ExprID::Mul, V, U>::getSizeAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Vector V, Scalar U>
    [[nodiscard]] auto operator*(const RSparseVector<V>& v, const U& x) noexcept {
        return SparseVectorExpr<ExprID::Mul, V, U>(v.getDerived(), x);
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator*(const U& x, const RSparseVector<V>& v) noexcept {
        return SparseVectorExpr<ExprID::Mul, V, U>(v * x);
    }
}

namespace Physica {
    template<ExprID ID, Vector V, Scalar U>
    class Traits<SparseVectorExpr<ID, V, U>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename V::ScalarType, U>::Type;
    };
}
