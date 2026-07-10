/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.h"

namespace Physica {
    template<class V> requires(std::remove_cvref_t<V>::isCompact())
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
        This& operator=(This&&) = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept { return 2; }
    };

    template<class V> requires(std::remove_cvref_t<V>::isCompact())
    auto ImagVector<V>::data_handle(this auto&& self) noexcept {
        return self.v.data()->imag_ptr();
    }
}
