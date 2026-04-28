/*
 * Copyright 2023-2026 Weibo He.
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

#include "../RValueVector.h"
#include "Physica/Core/Math/Algebra/Canonicalization.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/Utils/Unroller.h"

namespace Physica {
    template<Vector V1, Vector V2>
    class InnerDot {
        using This = InnerDot<V1, V2>;
        using T1 = V1::ScalarType;
        using T2 = V2::ScalarType;
        using T = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        using Tr = T::RealType;
    private:
        const V1& v1;
        const V2& v2;
    public:
        InnerDot(const V1& v1_, const V2& v2_);
        InnerDot(const This&) = delete;
        InnerDot(This&&) noexcept = delete;
        ~InnerDot() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        auto calc() const noexcept;
        T calc_mkl() const noexcept;
        CoDiff<T> calc_base() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    private:
        T calc_base_simd() const noexcept;
        T calc_base_simd_complex_real() const noexcept;
        T calc_base_fallback() const noexcept;
    };

    template<Vector V1, Vector V2>
    InnerDot<V1, V2>::InnerDot(const V1& v1_, const V2& v2_) : v1(v1_), v2(v2_) {
        assert(v1.getLength() == v2.getLength());
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc() const noexcept {
        if constexpr (T::isForwardDiff()) {
            if constexpr (!V1::isForwardDiff())
                return T(v2.values() * v1, v2.grads() * v1);
            else if constexpr (!V2::isForwardDiff())
                return T(v1.values() * v2, v1.grads() * v2);
            else {
                constexpr static int Order = T::Order - 1;
                return T(v1.values() * v2.values(), v2.template grads_mask<Order>() * v1.grads() + v1.template grads_mask<Order>() * v2.grads());
            }
        }
        else if constexpr (Internal::EnableMKL<V1, V2>::value)
            return calc_mkl();
        else
            return calc_base();
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc_base() const noexcept -> CoDiff<T> {
        if constexpr (T::isReverseDiff()) {
            auto& result = co_yield v1.values() * v2.values();
            if constexpr (ReverseDiff<V1>)
                v1.reverse(v2.values() * result.grad());
            if constexpr (ReverseDiff<V2>)
                v2.reverse(v1.values() * result.grad());
        }
        else if constexpr (V1::isFastPacket() && V2::isFastPacket()) {
            if constexpr (Internal::EnableSIMD<V1, V2>::value)
                co_return calc_base_simd();
            else if constexpr (T1::isComplex() && std::same_as<typename T1::RealType, T2>)
                co_return calc_base_simd_complex_real();
            else
                co_return calc_base_fallback();
        }
        else
            co_return calc_base_fallback();
    }

    template<Vector V1, Vector V2>
    __host__ __device__ consteval size_t InnerDot<V1, V2>::getSizeAtCompile() noexcept {
        return std::max(std::remove_cvref_t<V1>::getSizeAtCompile(), std::remove_cvref_t<V2>::getSizeAtCompile());
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc_base_simd() const noexcept -> T {
        using Pack = Internal::EnableSIMD<V1, V2>::PacketType;
        constexpr int Size = Pack::size();
        const size_t length = v1.getLength();
        const size_t to = length / Size * Size;
        size_t i = 0;
        auto buffer = Pack::zeros();
        auto it = zip(v1.view(), v2.view()).begin();
        for (; i < to; i += Size) {
            auto [lhs, rhs] = it + i;
            buffer = fma(lhs.template load<Size>(), rhs.template load<Size>(), buffer);
        }

        if (to != length) {
            auto [lhs, rhs] = it + i;
            const size_t count = length - i;
            buffer = fma(lhs.template load<Size>(count), rhs.template load<Size>(count), buffer);
        }
        return buffer.sum();
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc_base_simd_complex_real() const noexcept -> T {
        using Pack = BestPacket<T1, getSizeAtCompile()>::Type;
        using FullRealPack = Pack::FullRealType;
        const size_t length = v1.getLength();
        auto unroller = Unroller<Pack, HostDevAttr::NumUnrollDefault>();
        auto view1 = v1.view();
        auto view2 = v2.view();
        size_t i = unroller.loop_recursive([ite1 = view1.begin(), ite2 = view2.begin()](Pack buffer, size_t index) noexcept {
            auto p1 = (ite1 + index).template load<Pack::size()>().asReal();
            auto half = (ite2 + index).template load<FullRealPack::size() / 2>();
            auto p2 = FullRealPack(half, half).scatterRealImag();
            return Pack::asComplex(fma(p1, p2, buffer.asReal()));
        }, length);

        T result1 = unroller.sum().sum();
        T result2 = 0;
        for (; i < length; ++i) {
            auto x = v1.calc(i);
            auto y = v2.calc(i);
            result2.real() = fma(x.real(), y, result2.real());
            result2.imag() = fma(x.imag(), y, result2.imag());
        }
        return result1 + result2;
    }
    /**
     * Fallback if we do not know how to lower the inner dot
     */
    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc_base_fallback() const noexcept -> T {
        auto result = T(0);
        for(size_t i = 0; i < v1.getLength(); ++i)
            result += v1.calc(i) * v2.calc(i);
        return result;
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator*(const V1& v1, const V2& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        if constexpr (!canonicalized(v1, v2))
            return InnerDot<V2, V1>(v2, v1).calc();
        else
            return InnerDot<V1, V2>(v1, v2).calc();
    }
}

#ifdef PHYSICA_MKL
    #include "InnerDot_MKL.h"
#endif
