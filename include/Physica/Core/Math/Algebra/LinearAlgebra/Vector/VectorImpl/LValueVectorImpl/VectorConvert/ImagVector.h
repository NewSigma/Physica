/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica {
    template<class V> requires(std::remove_cvref_t<V>::isLValueVector() && !std::remove_cvref_t<V>::isStrided())
    class ImagVector<V> : public LValueVector<ImagVector<V>> {
        using This = ImagVector<V>;
        using Base = LValueVector<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<V> v;
    public:
        explicit ImagVector(V&& v_) : v(std::forward<V>(v_)) {}
        ImagVector(const This&) = default;
        ImagVector(This&&) noexcept = default;
        ~ImagVector() = default;
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
    auto ImagVector<V>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index out of range");
        return self.v[index].imag_ptr();
    }
}
