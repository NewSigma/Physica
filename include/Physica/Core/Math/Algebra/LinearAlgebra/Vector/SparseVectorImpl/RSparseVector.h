/*
 * Copyright 2024 WeiBo He.
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

namespace Physica::Core {
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

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RSparseVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        ResultType result(0);
        for (size_t i = 0; i < v1.getNumNonZero(); ++i) {
            const auto pair = v1.calcNonZero(i);
            result += ResultType(pair.second * v2.calc(pair.first));
        }
        return result;
    }

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RSparseVector<VectorType2>& v2) {
        return v2 * v1;
    }

    template<class Derived, class VectorType>
    void operator+=(LValueVector<Derived>& v1, const RSparseVector<VectorType>& v2) {
        using ResultType = typename Derived::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1[pair.first] += ResultType(pair.second);
        }
    }

    template<class Derived, class VectorType>
    void operator-=(LValueVector<Derived>& v1, const RSparseVector<VectorType>& v2) {
        using ResultType = typename Derived::ScalarType;
        for (size_t i = 0; i < v2.getNumNonZero(); ++i) {
            const auto pair = v2.calcNonZero(i);
            v1.getDerived()[pair.first] -= ResultType(pair.second);
        }
    }

    template<class Derived, class VectorType>
    inline void operator+=(ContinuousVector<Derived>& v1, const RSparseVector<VectorType>& v2) {
        static_cast<LValueVector<Derived>&>(v1) += v2;
    }

    template<class Derived, class VectorType>
    inline void operator-=(ContinuousVector<Derived>& v1, const RSparseVector<VectorType>& v2) {
        static_cast<LValueVector<Derived>&>(v1) -= v2;
    }
}

#include "SparseVectorExpr.h"
