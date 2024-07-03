/*
 * Copyright 2023-2024 WeiBo He.
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

#include <vectorclass/vectorclass.h>
#include <Physica/Utils/Container/Array/Array.h>
#include "Scalar.h"
#include "SIMDImpl/Instruset.h"

namespace Physica::Core {
    template<class T, size_t Size> class SIMD;
    template<class T, size_t Size> class BoolSIMD;
    template<class ScalarType, unsigned int Order> class DiffTracer;

    namespace Internal {
        /**
         * Find the best packet for a linear storage
         */
        template<class ScalarType, size_t Length>
        class BestPacket {
            using PlainScalar = typename ScalarType::PlainScalar;
            constexpr static bool isSinglePrec = ScalarType::Option == Float;
            constexpr static bool isComplex = ScalarType::isComplex;
            constexpr static bool isForward = ScalarType::isForwardDiff;
            constexpr static bool isDynamic = Length == Dynamic;
            constexpr static size_t size128 = (isSinglePrec ? 4 : 2) / (isComplex ? 2 : 1) / (isForward ? 2 : 1);
            constexpr static size_t size256 = (isSinglePrec ? 8 : 4) / (isComplex ? 2 : 1) / (isForward ? 2 : 1);
            constexpr static size_t size512 = (isSinglePrec ? 16 : 8) / (isComplex ? 2 : 1) / (isForward ? 2 : 1);
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
            constexpr static size_t Size = (isComplex || isForward) ? 1 : (isDynamic ? BiggestSize : Size3);
            using Type = typename std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
        };
    }

    template<class ScalarType, size_t Size>
    class SIMD : private Traits<SIMD<ScalarType, Size>>::BaseType {
        using This = SIMD<ScalarType, Size>;
        using BoolSIMDType = typename Traits<This>::BoolSIMDType;
        using PlainScalar = typename ScalarType::PlainScalar;
        using TrivialType = typename PlainScalar::TrivialType;
        constexpr static bool isForward = ScalarType::isForwardDiff;
    public:
        using Base = typename Traits<This>::BaseType;
        using PlainPacket = This;
    public:
        SIMD() = default;
        explicit SIMD(ScalarType s) : Base(s.getTrivial()) {}
        explicit SIMD(Base value) : Base(value) {}
        using Base::Base;
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
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
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] Base& getImpl() noexcept { return *this; }
        [[nodiscard]] const Base& getImpl() const noexcept { return *this; }
    };

    template<class PlainScalar, size_t Size>
    class SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size> : public SIMD<PlainScalar, Size> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Invalid template param");
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse, 1>;
        using This = SIMD<ScalarType, Size>;
        using Base = SIMD<PlainScalar, Size>;
        using BoolSIMDType = typename Traits<This>::BoolSIMDType;
        using TrivialType = typename PlainScalar::TrivialType;
    public:
        using PlainPacket = Base;
        using TracerType = DiffTracer<PlainScalar, 1>;
    private:
        ScalarType headNode;
    public:
        SIMD() = default;
        explicit SIMD(int i) : SIMD(PlainScalar(i)) {}
        explicit SIMD(PlainScalar s);
        explicit SIMD(ScalarType s);
        SIMD(Base base, ScalarType headNode_);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] inline ScalarType operator[](int i) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        using Base::operator>;
        using Base::operator<;
        using Base::operator>=;
        using Base::operator<=;
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(int n, const ScalarType* p);
        inline void store(ScalarType* p) const;
        inline void store_partial(int n, ScalarType* p) const;
        [[nodiscard]] inline ScalarType horizontal_add() const;
        inline void swap(SIMD& __restrict obj) noexcept;
        /* Getters */
        using Base::size;
        using Base::getImpl;
        [[nodiscard]] ScalarType getHeadNode() const noexcept { return headNode; }
        [[nodiscard]] PlainScalar* value_ptr() const noexcept { return headNode.value_ptr(); }
        [[nodiscard]] PlainScalar* grad_ptr() const noexcept { return headNode.grad_ptr(); }
    private:
        using Base::insert; //Insert a scalar may lead to incontineous memory, which harms performance
        [[nodiscard]] static bool checkContinuous(int n, const ScalarType* p);
    };

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> operator*(const SIMD<ScalarType, Size>& packet, const ScalarType& scalar) {
        return SIMD<ScalarType, Size>(packet.getImpl() * scalar.getValue().getTrivial());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> operator*(const ScalarType& scalar, const SIMD<ScalarType, Size>& packet) {
        return packet * scalar;
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> mul_add(
            const SIMD<ScalarType, Size>& a,
            const SIMD<ScalarType, Size>& b,
            const SIMD<ScalarType, Size>& c);

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> nmul_add(
            const SIMD<ScalarType, Size> a,
            const SIMD<ScalarType, Size> b,
            const SIMD<ScalarType, Size> c);

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> mul_sub(
            const SIMD<ScalarType, Size> a,
            const SIMD<ScalarType, Size> b,
            const SIMD<ScalarType, Size> c);
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Size>
    class Traits<SIMD<T, Size>> {
    public:
        using ScalarType = T;

        static_assert(ScalarType::Option == Float || ScalarType::Option == Double, "[Error]: Unsupported float type");
        static_assert(!ScalarType::isComplex, "[Error]: The main template targets on real scalar");
        static_assert(Size % 2 == 0 && Size <= 16, "[Error]: Invalid Size");
    private:
        using PlainScalar = typename ScalarType::PlainScalar;
        constexpr static bool isForward = ScalarType::isForwardDiff;
        constexpr static bool isSinglePrec = ScalarType::Option == Float;
        constexpr static size_t PlainSize = isForward ? (Size * 2) : Size;
        static_assert(!isForward, "[Error]: SIMD for forward autodiff is not implemented");
            
        using Size2Type = typename std::conditional<isSinglePrec, void, Vec2d>::type;
        using Size4Type = typename std::conditional<isSinglePrec, Vec4f, Vec4d>::type;
        using Size8Type = typename std::conditional<isSinglePrec, Vec8f, Vec8d>::type;
        using Size16Type = typename std::conditional<isSinglePrec, Vec16f, void>::type;
        using Type1 = typename std::conditional<PlainSize == 2, Size2Type, Size4Type>::type;
        using Type2 = typename std::conditional<PlainSize == 8, Size8Type, Size16Type>::type;
    public:
        using BaseType = typename std::conditional<PlainSize <= 4, Type1, Type2>::type;
        using BoolSIMDType = BoolSIMD<ScalarType, Size>;
    };
}

#include "SIMDImpl/BoolSIMD.h"
#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/ElementaryFunction.h"
