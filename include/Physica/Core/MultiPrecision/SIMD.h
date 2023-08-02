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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
    #include "vectorclass/vectorclass.h"
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#include "Scalar.h"
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    template<class T, size_t Size> class SIMD;
    template<class T, size_t Size> class BoolSIMD;

    namespace Internal {
        class Instrset {
        public:
            [[nodiscard]] constexpr static bool hasAVX() {
            #ifdef PHYSICA_AVX
                return true;
            #else
                return false;
            #endif
            }

            [[nodiscard]] constexpr static bool hasAVX2() {
            #ifdef PHYSICA_AVX2
                return true;
            #else
                return false;
            #endif
            }

            [[nodiscard]] constexpr static bool hasAVX512() {
            #ifdef PHYSICA_AVX512
                return true;
            #else
                return false;
            #endif
            }
        };

        template<class T, size_t Size>
        class Traits<SIMD<T, Size>> {
        public:
            using ScalarType = T;
        private:
            static_assert(ScalarType::option == Float || ScalarType::option == Double, "[Error]: Unsupported float type");
            static_assert(!ScalarType::isComplex, "[Error]: The main template targets on real scalar");
            static_assert(!ScalarType::isDifferentiable, "[Error]: The main template is not differentiable");
            static_assert(Size % 2 == 0 && Size <= 16, "[Error]: Invalid Size");
            constexpr static bool isSinglePrec = ScalarType::option == Float;
            
            using Size2Type = typename std::conditional<isSinglePrec, void, Vec2d>::type;
            using Size4Type = typename std::conditional<isSinglePrec, Vec4f, Vec4d>::type;
            using Size8Type = typename std::conditional<isSinglePrec, Vec8f, Vec8d>::type;
            using Size16Type = typename std::conditional<isSinglePrec, Vec16f, void>::type;
            using Type1 = typename std::conditional<Size == 2, Size2Type, Size4Type>::type;
            using Type2 = typename std::conditional<Size == 8, Size8Type, Size16Type>::type;
        public:
            using BaseType = typename std::conditional<Size <= 4, Type1, Type2>::type;
            using BoolSIMDType = BoolSIMD<ScalarType, Size>;
        };
        /**
         * Find the best packet for a linear storage
         */
        template<class ScalarType, size_t Length>
        class BestPacket {
            constexpr static bool isSinglePrec = ScalarType::option == Float;
            constexpr static bool isComplex = ScalarType::isComplex;
            constexpr static bool isDifferentiable = ScalarType::isDifferentiable;
            constexpr static bool isDynamic = Length == Utils::Dynamic;
            constexpr static size_t size128 = (isSinglePrec ? 4 : 2) / (isComplex ? 2 : 1) / (isDifferentiable ? 2 : 1);
            constexpr static size_t size256 = (isSinglePrec ? 8 : 4) / (isComplex ? 2 : 1) / (isDifferentiable ? 2 : 1);
            constexpr static size_t size512 = (isSinglePrec ? 16 : 8) / (isComplex ? 2 : 1) / (isDifferentiable ? 2 : 1);
            constexpr static bool support128 = INSTRSET >= 2;
            constexpr static bool support256 = INSTRSET >= 7 && support128;
            constexpr static bool support512 = INSTRSET >= 9 && support256;
            constexpr static bool use128 = support128 && Length >= size128 && size128 != 0;
            constexpr static bool use256 = support256 && Length >= size256 && size256 != 0;
            constexpr static bool use512 = support512 && Length >= size512 && size512 != 0;
            constexpr static size_t Size1 = use128 ? size128 : 1;
            constexpr static size_t Size2 = use256 ? size256 : Size1;
            constexpr static size_t Size3 = use512 ? size512 : Size2;
            constexpr static size_t BiggestSize1 = support128 ? size128 : 1;
            constexpr static size_t BiggestSize2 = support256 ? size256 : BiggestSize1;
            constexpr static size_t BiggestSize = support512 ? size512 : BiggestSize2;
        public:
            constexpr static size_t Size = isDynamic ? BiggestSize : Size3;
            using Type = typename std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
        };
    }

    template<class ScalarType, size_t Size>
    class SIMD : private Internal::Traits<SIMD<ScalarType, Size>>::BaseType {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid template param");
        using This = SIMD<ScalarType, Size>;
        using Traits = Internal::Traits<This>;
        using BoolSIMDType = typename Traits::BoolSIMDType;
        using TrivialType = typename ScalarType::TrivialType;
    public:
        using Base = typename Traits::BaseType;
    public:
        SIMD() = default;
        explicit SIMD(Base value) : Base(value) {}
        using Base::Base;
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] inline ScalarType operator[](int index) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        [[nodiscard]] inline BoolSIMDType operator>(const SIMD& other) const;
        [[nodiscard]] inline BoolSIMDType operator<(const SIMD& other) const;
        [[nodiscard]] inline BoolSIMDType operator>=(const SIMD& other) const { return !(*this < other); }
        [[nodiscard]] inline BoolSIMDType operator<=(const SIMD& other) const { return !(*this > other); }
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(int n, const ScalarType* p);
        inline void store(ScalarType* p) const;
        inline void store_partial(int n, ScalarType* p) const;
        inline void insert(int index, const ScalarType& value);
        [[nodiscard]] inline ScalarType horizontal_add() const;
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] Base& getImpl() noexcept { return *this; }
        [[nodiscard]] const Base& getImpl() const noexcept { return *this; }
    };

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> operator*(const SIMD<ScalarType, Size>& packet, const ScalarType& scalar) {
        return SIMD<ScalarType, Size>(packet.getImpl() * scalar.getTrivial());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> operator*(const ScalarType& scalar, const SIMD<ScalarType, Size>& packet) {
        return packet * scalar;
    }
}

#include "SIMDImpl/BoolSIMD.h"
#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/ElementaryFunction.h"
