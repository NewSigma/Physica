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
        static_assert(V::isComplex, "[Error]: Unnecessary conjugate call on real vector");

        using This = Conjugate<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
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
        void assign(Vector auto&& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Vector auto& v) const noexcept;

        [[nodiscard]] T calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] Tv calc_value(size_t index) const { return vec.calc_value(index).conjugate(); }
        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    template<ExecutePolicy P>
    void Conjugate<V>::assign(Vector auto&& v) const noexcept {
        assign_base<P>(v);
    }

    template<Vector V>
    template<ExecutePolicy P>
    void Conjugate<V>::assign_base(Vector auto& v) const noexcept {
        Base::template assign<P>(v);
    }

    template<Vector V>
    template<Packet Pack>
    Pack Conjugate<V>::packet(size_t index) const noexcept {
        return vec.template packet<Pack>(index).conjugate();
    }

    template<Vector V>
    template<Packet Pack>
    Pack Conjugate<V>::packetPartial(size_t index, size_t count) const noexcept {
        return vec.template packetPartial<Pack>(index, count).conjugate();
    }
}

namespace Physica {
    template<Vector V>
    class Traits<Conjugate<V>> : public Traits<V> {};
}

#ifdef PHYSICA_MKL
    #include "Conjugate_MKL.h"
#endif
