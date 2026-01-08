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
        [[nodiscard]] size_t getNumNonZero() const noexcept { return Base::getDerived().getNumNonZero(); }
    };

    template<Vector V1, Vector V2>
    Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type
    operator*(const RSparseVector<V1>& v1, const V2& v2) {
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
    Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type
    operator*(const V1& v1, const RSparseVector<V2>& v2) {
        return v2 * v1;
    }

    template<Vector V, Vector U>
    void operator+=(V& v1, const U& v2) requires Sparse<U> {
        using ResultType = V::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1[pair.first] += ResultType(pair.second);
        }
    }

    template<Vector V, Vector U>
    void operator-=(V& v1, const U& v2) requires Sparse<U> {
        using ResultType = V::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1[pair.first] -= ResultType(pair.second);
        }
    }
}

#include "SparseVectorExpr.h"
