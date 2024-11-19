/*
 * Copyright 2023-2024 Weibo He.
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
    template<Vector T>
    class ReverseVector<T> final : public RValueVector<ReverseVector<T>> {
        using This = ReverseVector<T>;
        using Base = RValueVector<This>;

        T& v;
    public:
        using typename Base::ScalarType;
    public:
        explicit ReverseVector(T& v_) : v(v_) {}
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
    template<Vector T>
    class Traits<Core::ReverseVector<T>> : public Traits<T> {};
}
