/*
 * Copyright 2020-2025 Weibo He.
 *
 * This file is part of Physica.

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

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <concepts>
#include <coroutine>
#include "Physica/Macro.h"
#include "Physica/Core/Utils/Type.h"

namespace Physica {
    namespace Internal {
        template<class T>
        struct IsScalarRef {
            constexpr static bool value = false;
        };

        template<class T>
        struct IsScalarRef<ScalarRef<T>> {
            constexpr static bool value = true;
        };
    }

    enum FloatPrec : char {
        Float16 = 0,
        Float32 = 1,
        Float64 = 2,
        FloatMP = 3,

        Half = Float16,
        Float = Float32,
        Double = Float64
    };

    enum class DiffMode : char {
        Forward,
        Reverse
    };

    // MP = MultiPrecision
    using MPUnit = std::conditional<PhysicaWordSize == 64, uint64_t, uint32_t>::type;
    using SignedMPUnit = std::conditional<PhysicaWordSize == 64, int64_t, int32_t>::type;
    constexpr static unsigned int MPUnitWidth = PhysicaWordSize;
    constexpr static MPUnit MPUnitMax = std::numeric_limits<MPUnit>::max();
    constexpr static MPUnit MPUnitHighestBitMask = static_cast<MPUnit>(1) << (MPUnitWidth - 1);
    constexpr static MPUnit MPUnitLowerMask = MPUnitMax >> (MPUnitWidth / 2);

    template<class ScalarType>
    class ScalarPtr;
    template<class Derived>
    class ScalarBase;
    template<class Derived>
    class SIMDBase;
    /**
     * \class Real is a advanced float type that supports multiple precision
     */
    template<FloatPrec Prec = Float64>
    class Real;
    template<class T>
    class Complex;
    template<class ScalarType, DiffMode Mode, int Order = 1>
    class Diff;

    template<class T>
    concept Scalar = std::derived_from<std::remove_cvref_t<T>, ScalarBase<std::remove_cvref_t<T>>>
                  || std::derived_from<std::remove_cvref_t<T>, typename std::remove_cvref_t<T>::ScalarType>
                  || Internal::IsScalarRef<std::remove_cvref_t<T>>::value;

    template<class T>
    concept Packet = Scalar<T> || std::derived_from<std::remove_cvref_t<T>, SIMDBase<std::remove_cvref_t<T>>>;

    template<class T>
    concept ForwardDiff = std::remove_cvref_t<T>::ScalarType::isForwardDiff;

    template<class T>
    concept ReverseDiff = std::remove_cvref_t<T>::ScalarType::isReverseDiff;

    template<class T>
    concept Diffable = std::remove_cvref_t<T>::ScalarType::isDiffable;

    template<class T>
    class DiffCoro;

    template<class T>
    class IsCoDiff {
        template<class U>
        struct Impl {
            constexpr static bool value = false;
        };

        template<class U>
        struct Impl<DiffCoro<U>> {
            constexpr static bool value = true;
        };
    public:
        constexpr static bool value = Impl<std::remove_cvref_t<T>>::value;
    };

    template<class T>
    struct remove_codiff {
        using Type = T;
    };

    template<class T>
    using remove_codiff_t = remove_codiff<T>::Type;

    template<class T>
    struct remove_codiff<DiffCoro<T>> {
        using Type = T;
    };

    template<class T>
    using CoDiff = std::conditional<std::is_void<T>::value || ReverseDiff<T>
                 , DiffCoro<remove_codiff_t<typename remove_cvref<T>::Type>>
                 , typename remove_cvref<T>::Type>::type;

    template<class T> requires(std::is_reference<T>::value)
    using LazyDestroy = std::conditional<std::is_rvalue_reference<T>::value, std::remove_reference_t<T>, std::add_lvalue_reference_t<T>>::type;

    namespace Internal {
        /**
         * \class BinaryScalarOpRtnTy returns the minimal type that can exactly represent the two input scalars.
         */
        template<Scalar T1, Scalar T2>
        class BinaryScalarOpRtnTy {
            static_assert(std::is_class<T1>::value);
            static_assert(std::is_class<T2>::value);

            constexpr static FloatPrec Prec = std::max(T1::Prec, T2::Prec);
            constexpr static bool isComplex = T1::isComplex || T2::isComplex;
            constexpr static bool isDiffable1 = T1::isDiffable;
            constexpr static bool isDiffable2 = T2::isDiffable;
            constexpr static bool isDiffable = isDiffable1 || isDiffable2;

            constexpr static DiffMode Mode = T1::isDiffable ? T1::Mode : T2::Mode;
            constexpr static DiffMode Mode1 = T2::isDiffable ? T2::Mode : T1::Mode;
            static_assert(Mode == Mode1, "[Error]: Operation between differentiable scalars with different mode is not well defined");

            constexpr static int Order1 = T1::Order;
            constexpr static int Order2 = T2::Order;
            constexpr static bool UseMixOrder = isDiffable1 && isDiffable2 && (Order1 != Order2);
            constexpr static int Order = UseMixOrder ? std::min(Order1, Order2) : std::max(Order1, Order2);
            static_assert(Mode != DiffMode::Reverse || !UseMixOrder, "[Error]: Reverse mode does not support mixed order");

            using Type0 = Real<Prec>;
            using Type1 = std::conditional<isComplex, Complex<Type0>, Type0>::type;
            using Type2 = std::conditional<isDiffable, Diff<Type1, Mode, Order>, Type1>::type;
        public:
            using Type = Type2;
        };
    }

    template<>
    class DiffCoro<void> {
        using This = DiffCoro<void>;
        class Promise {
        public:
            auto get_return_object() { return std::coroutine_handle<Promise>::from_promise(*this); };
            static std::suspend_never initial_suspend() noexcept { return {}; }
            static std::suspend_always final_suspend() noexcept { return {}; }
            void return_void() noexcept {}
            [[noreturn]] static void unhandled_exception() { throw; }
        };
    public:
        using promise_type = Promise;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        DiffCoro() = default;
        DiffCoro(std::coroutine_handle<Promise> handle_) noexcept : handle(handle_) {}
        DiffCoro(const This&) = delete;
        DiffCoro(This&& other) noexcept : handle(other.handle) { other.handle = nullptr; }
        ~DiffCoro() {
            if (!handle.done()) {
                handle.resume();
                assert(handle.done() && "[Error]: Invalid reverse diff");
                handle.destroy();
                handle = nullptr;
            }
        }
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            std::swap(handle, obj.handle);
        }
    };
}
