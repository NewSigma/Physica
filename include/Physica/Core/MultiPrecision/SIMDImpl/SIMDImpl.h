/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline ScalarType SIMD<ScalarType, Size>::operator[](int index) const {
        if constexpr (isForward)
            return ScalarType(Base::operator[](index * 2), Base::operator[](index * 2 + 1));
        else
            return ScalarType(Base::operator[](index));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator+(const SIMD& other) const {
        return SIMD(getImpl() + other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator-(const SIMD& other) const {
        return SIMD(getImpl() - other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator*(const SIMD& other) const {
        return SIMD(getImpl() * other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator/(const SIMD& other) const {
        return SIMD(getImpl() / other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator-() const {
        return SIMD(-getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline typename SIMD<ScalarType, Size>::BoolSIMDType
    SIMD<ScalarType, Size>::operator>(const SIMD& other) const {
        return BoolSIMDType(getImpl() > other.getImpl());
    }
    
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline typename SIMD<ScalarType, Size>::BoolSIMDType
    SIMD<ScalarType, Size>::operator<(const SIMD& other) const {
        return BoolSIMDType(getImpl() < other.getImpl());
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::load(const ScalarType* p) {
        Base::load(reinterpret_cast<const TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::load_partial(int n, const ScalarType* p) {
        Base::load_partial(n, reinterpret_cast<const TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::store(ScalarType* p) const {
        Base::store(reinterpret_cast<TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::store_partial(int n, ScalarType* p) const {
        Base::store_partial(n, reinterpret_cast<TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::insert(int index, const ScalarType& value) {
        if constexpr (isForward) {
            Base::insert(index * 2, value.getValue().getTrivial());
            Base::insert(index * 2 + 1, value.getTangent().getTrivial());
        }
        else
            Base::insert(index, value.getTrivial());
    }

    template<class ScalarType, size_t Size>
    inline ScalarType SIMD<ScalarType, Size>::horizontal_add() const {
        using Physica::horizontal_add;
        return horizontal_add(getImpl());
    }
}
