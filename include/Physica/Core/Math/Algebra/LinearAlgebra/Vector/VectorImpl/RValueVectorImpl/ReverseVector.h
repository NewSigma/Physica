/*
 * Copyright 2023-2024 WeiBo He.
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

namespace Physica::Core {
    template<class VectorType>
    class ReverseVector final : public RValueVector<ReverseVector<VectorType>> {
        using This = ReverseVector<VectorType>;
        using Base = RValueVector<This>;

        VectorType& v;
    public:
        using typename Base::ScalarType;
    public:
        explicit ReverseVector(RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        ReverseVector(const ReverseVector&) = delete;
        ReverseVector(ReverseVector&&) noexcept = delete;
        ~ReverseVector() = default;
        /* Operators */
        ReverseVector& operator=(const ReverseVector&) = delete;
        ReverseVector& operator=(ReverseVector&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return v.calc(getLength() - index - 1); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };
}

namespace Physica {
    template<class VectorType>
    class Traits<Core::ReverseVector<VectorType>> : public Traits<VectorType> {};
}
