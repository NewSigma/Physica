/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
    template<class V>
    class RealVectorR : public RValueVector<RealVectorR<V>> {
        using This = RealVectorR<V>;
        using Base = RValueVector<This>;
        static_assert(V::isComplex, "[Error]: Unnecessary real() call on real vector");
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        explicit RealVectorR(const V& v_) : v(v_) {}
        RealVectorR(const This&) = default;
        RealVectorR(This&&) noexcept = default;
        ~RealVectorR() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).real(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s).real(); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };
}

namespace Physica {
    template<class V>
    class Traits<RealVectorR<V>> {
    public:
        using ScalarType = V::ScalarType::RealType;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = Traits<V>::FastPacket;
    };
}
