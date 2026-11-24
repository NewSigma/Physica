/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
    template<class Derived>
    class RSparseVector : public RValueVector<Derived> {
        using Base = RValueVector<Derived>;
    public:
        using typename Base::ScalarType;
        using NonZeroPair = std::pair<size_t, ScalarType>;
    public:
        /* Getters */
        [[nodiscard]] NonZeroPair calcNonZero(size_t index) const { return Base::getDerived().calcNonZero(index); }
        [[nodiscard]] constexpr size_t getNumNonZero() const noexcept { return Base::getDerived().getNumNonZero(); }
    };

    template<Vector V1, Vector V2>
    auto operator*(const RSparseVector<V1>& v1, const V2& v2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        ResultType result(0);
        for (size_t i = 0; i < v1.getNumNonZero(); ++i) {
            const auto pair = v1.calcNonZero(i);
            result += ResultType(pair.second * v2.calc(pair.first));
        }
        return result;
    }

    template<Vector V1, Vector V2>
    void operator+=(V1& v1, const V2& v2) requires(V2::isSparse()) {
        assert(v1.getLength() == v2.getLength());
        using ResultType = V1::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1[pair.first] += ResultType(pair.second);
        }
    }

    template<Vector V1, Vector V2>
    void operator-=(V1& v1, const V2& v2) requires(V2::isSparse()) {
        assert(v1.getLength() == v2.getLength());
        using ResultType = V1::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1[pair.first] -= ResultType(pair.second);
        }
    }
}

#include "SparseVectorExpr.h"
