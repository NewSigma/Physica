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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica {
    template<class V>
    class ImagVectorL : public LValueVector<ImagVectorL<V>> {
        using This = ImagVectorL<V>;
        using Base = LValueVector<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<V> v;
    public:
        explicit ImagVectorL(V&& v_) : v(std::forward<V>(v_)) {}
        ImagVectorL(const This&) = default;
        ImagVectorL(This&&) noexcept = default;
        ~ImagVectorL() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<class V>
    auto ImagVectorL<V>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index out of range");
        return self.v[index].imag_ptr();
    }
}

namespace Physica {
    template<class V>
    class Traits<ImagVectorL<V>> : public Traits<ImagVectorR<V>> {};
}
