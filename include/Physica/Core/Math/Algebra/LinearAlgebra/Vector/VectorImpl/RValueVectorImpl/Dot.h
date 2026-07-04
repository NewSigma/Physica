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
    template<Vector LHS, Vector RHS>
    class Dot {
        using This = Dot<LHS, RHS>;
        using T1 = std::remove_cvref_t<LHS>::ScalarType;
        using T2 = std::remove_cvref_t<RHS>::ScalarType;
        using T = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        using Tr = T::RealType;
    private:
        decay_rvalue_t<LHS> lhs;
        decay_rvalue_t<RHS> rhs;
    public:
        Dot(LHS lhs_, RHS rhs_);
        Dot(const This&) = default;
        Dot(This&&) noexcept = default;
        ~Dot() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] auto calc(this auto&&) noexcept;
        [[nodiscard]] T calc_mkl() const noexcept;
        [[nodiscard]] CoDiff<T> calc_base(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    private:
        [[nodiscard]] T calc_base_simd() const noexcept;
        [[nodiscard]] T calc_base_simd_complex_real() const noexcept;
        [[nodiscard]] T calc_base_fallback() const noexcept;

        [[nodiscard]] constexpr bool useMKL() const noexcept;
    };

    template<Vector LHS, Vector RHS>
    Dot<LHS, RHS>::Dot(LHS lhs_, RHS rhs_) : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        assert(lhs.getLength() == rhs.getLength());
    }

    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (T::isForwardDiff()) {
            const auto& lhs = self.lhs;
            const auto& rhs = self.rhs;
            if constexpr (!lhs.isForwardDiff())
                return T(rhs.values() * lhs, rhs.grads() * lhs);
            else if constexpr (!rhs.isForwardDiff())
                return T(lhs.values() * rhs, lhs.grads() * rhs);
            else {
                constexpr static int Order = T::Order - 1;
                return T(lhs.values() * rhs.values(), rhs.template grads_mask<Order>() * lhs.grads() + lhs.template grads_mask<Order>() * rhs.grads());
            }
        }
        else if constexpr (HasMKL() && Internal::EnableLAPACK<LHS, RHS>::value)
            return self.useMKL() ? self.calc_mkl() : std::forward<Self>(self).calc_base();
        else
            return std::forward<Self>(self).calc_base();
    }

    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc_base(this auto&& self) noexcept -> CoDiff<T> {
        using Self = decltype(self);
        if constexpr (T::isReverseDiff()) {
            decltype(auto) dot = decay_rvalue(std::forward<Self>(self));
            const auto& lhs = dot.lhs;
            const auto& rhs = dot.rhs;
            auto& result = co_yield lhs.values() * rhs.values();
            if constexpr (ReverseDiff<LHS>)
                lhs.reverse(rhs.values() * result.grad());
            if constexpr (ReverseDiff<RHS>)
                rhs.reverse(lhs.values() * result.grad());
        }
        else if constexpr (self.lhs.isFastPacket() && self.rhs.isFastPacket()) {
            if constexpr (Internal::EnableSIMD<LHS, RHS>::value)
                co_return self.calc_base_simd();
            else if constexpr (T1::isComplex() && std::same_as<typename T1::RealType, T2>)
                co_return self.calc_base_simd_complex_real();
            else
                co_return self.calc_base_fallback();
        }
        else
            co_return self.calc_base_fallback();
    }

    template<Vector LHS, Vector RHS>
    __host__ __device__ consteval size_t Dot<LHS, RHS>::getSizeAtCompile() noexcept {
        return std::max(std::remove_cvref_t<LHS>::getSizeAtCompile(), std::remove_cvref_t<RHS>::getSizeAtCompile());
    }

    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc_base_simd() const noexcept -> T {
        using Pack = Internal::EnableSIMD<LHS, RHS>::PacketType;
        constexpr int Size = Pack::size();
        const size_t length = lhs.getLength();
        const size_t to = length / Size * Size;
        size_t i = 0;
        auto buffer = Pack::zeros();
        auto it = zip(lhs.view(), rhs.view()).begin();
        for (; i < to; i += Size) {
            auto [v1, v2] = it + i;
            buffer = fma(v1.template load<Size>(), v2.template load<Size>(), buffer);
        }

        if (to != length) {
            auto [v1, v2] = it + i;
            const size_t count = length - i;
            buffer = fma(v1.template load<Size>(count), v2.template load<Size>(count), buffer);
        }
        return buffer.sum();
    }

    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc_base_simd_complex_real() const noexcept -> T {
        using Pack = BestPacket<T1, getSizeAtCompile()>::Type;
        using FullRealPack = Pack::FullRealType;
        const size_t length = lhs.getLength();
        auto unroller = Unroller<Pack, HostDevAttr::NumUnrollDefault>();
        auto view1 = lhs.view();
        auto view2 = rhs.view();
        size_t i = unroller.loop_recursive([ite1 = view1.begin(), ite2 = view2.begin()](Pack buffer, size_t index) noexcept {
            auto p1 = (ite1 + index).template load<Pack::size()>().asReal();
            auto half = (ite2 + index).template load<FullRealPack::size() / 2>();
            auto p2 = FullRealPack(half, half).scatterRealImag();
            return Pack::asComplex(fma(p1, p2, buffer.asReal()));
        }, length);

        T result1 = unroller.sum().sum();
        T result2 = 0;
        for (; i < length; ++i) {
            auto x = lhs.calc(i);
            auto y = rhs.calc(i);
            result2.real() = fma(x.real(), y, result2.real());
            result2.imag() = fma(x.imag(), y, result2.imag());
        }
        return result1 + result2;
    }
    /**
     * Fallback if we do not know how to lower the dot
     */
    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc_base_fallback() const noexcept -> T {
        auto result = T(0);
        for(size_t i = 0; i < lhs.getLength(); ++i)
            result += lhs.calc(i) * rhs.calc(i);
        return result;
    }

    template<Vector LHS, Vector RHS>
    constexpr bool Dot<LHS, RHS>::useMKL() const noexcept {
        constexpr int Threshold = HostDevAttr::ThresholdDot_MKL;
        if constexpr (getSizeAtCompile() == Dynamic)
            return Threshold < lhs.getLength();
        else
            return Threshold < getSizeAtCompile();
    }
    /**
     * \returns a Dot object with proper canonicalization, while operator* is syntactic suger for it
     */
    template<Vector LHS, Vector RHS>
    [[nodiscard]] auto dot(LHS&& lhs, RHS&& rhs) noexcept requires(!DeviceObj<LHS> && !DeviceObj<RHS>) {
        if constexpr (!canonicalized(lhs, rhs))
            return Dot<RHS&&, LHS&&>(std::forward<RHS>(rhs), std::forward<LHS>(lhs));
        else
            return Dot<LHS&&, RHS&&>(std::forward<LHS>(lhs), std::forward<RHS>(rhs));
    }

    [[nodiscard]] __host__ __device__ auto operator*(Vector auto&& lhs, Vector auto&& rhs) noexcept {
        return dot(std::forward<decltype(lhs)>(lhs), std::forward<decltype(rhs)>(rhs)).calc();
    }
}

#ifdef PHYSICA_MKL
    #include "Dot_MKL.h"
#endif
