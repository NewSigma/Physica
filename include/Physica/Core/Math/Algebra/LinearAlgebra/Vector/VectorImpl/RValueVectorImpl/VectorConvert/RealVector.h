/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
    template<class V>
    class RealVector : public RValueVector<RealVector<V>> {
        using This = RealVector<V>;
        using Base = RValueVector<This>;
        static_assert(std::remove_cvref_t<V>::isComplex(), "[Error]: Unnecessary real() call on real vector");
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<V> v;
    public:
        explicit RealVector(V&& v_) : v(std::forward<V>(v_)) {}
        RealVector(const This&) = default;
        RealVector(This&&) noexcept = default;
        ~RealVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).real(); }

        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<class V>
    decltype(auto) RealVector<V>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v).values().reals();
    }

    template<class V>
    __host__ __device__ consteval size_t RealVector<V>::getSizeAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }
}

namespace Physica {
    template<class V>
    class Traits<RealVector<V>> {
        using V1 = std::remove_cvref_t<V>;
    public:
        using ScalarType = V1::ScalarType::RealType;
    };
}
