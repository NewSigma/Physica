/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica {
    template<class V> requires(std::remove_cvref_t<V>::isLValueVector() && !std::remove_cvref_t<V>::isStrided())
    class RealVector<V> : public LValueVector<RealVector<V>> {
        using This = RealVector<V>;
        using Base = LValueVector<This>;
    protected:
        using typename Base::T;
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
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<class V> requires(std::remove_cvref_t<V>::isLValueVector() && !std::remove_cvref_t<V>::isStrided())
    auto RealVector<V>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index out of range");
        return self.v[index].real_ptr();
    }
}
