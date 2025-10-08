/*
 * Copyright 2022-2025 Weibo He.
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

namespace Physica {
    template<Vector V>
    class Conjugate<V> : public RValueVector<Conjugate<V>> {
        using This = Conjugate<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& vec;
    public:
        explicit Conjugate(const V& vec_) : vec(vec_) {}
        Conjugate(const This&) = delete;
        Conjugate(This&&) = delete;
        ~Conjugate() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        /* Getters */
        [[nodiscard]] T calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] Tv calc_value(size_t index) const { return vec.calc_value(index).conjugate(); }
        [[nodiscard]] size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    template<ExecutePolicy P>
    void Conjugate<V>::assign(Vector auto& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target[i] = calc(i);
    }
}

namespace Physica {
    template<Vector V>
    class Traits<Conjugate<V>> : public Traits<V> {};
}
