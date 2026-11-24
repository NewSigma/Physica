/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.h"

namespace Physica {
    template<class V> requires(std::remove_cvref_t<V>::isStrided())
    class ImagVector<V> : public StridedVector<ImagVector<V>> {
        using This = ImagVector<V>;
        using Base = StridedVector<This>;
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
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] constexpr size_t getStride() const noexcept { return 2 * v.getStride(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept { return 2 * std::remove_cvref_t<V>::getStrideAtCompile(); }
    };

    template<class V> requires(std::remove_cvref_t<V>::isStrided())
    auto ImagVector<V>::data_handle(this auto&& self) noexcept {
        return self.v[0].imag_ptr();
    }
}
