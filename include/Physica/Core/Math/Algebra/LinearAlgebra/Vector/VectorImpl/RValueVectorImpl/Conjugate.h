/*
 * Copyright 2022-2026 Weibo He.
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
        LazyDestroy<V> vec;
    public:
        explicit Conjugate(V&& vec_) : vec(std::forward<V>(vec_)) {}
        Conjugate(const This&) = default;
        Conjugate(This&&) = default;
        ~Conjugate() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] T calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] Tv calc_value(size_t index) const { return vec.calc_value(index).conjugate(); }
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;

        [[nodiscard]] decltype(auto) conjugate(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return vec.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Vector V>
    template<ExecutePolicy P>
    void Conjugate<V>::assign(Vector auto&& v) const noexcept {
        Base::template assign_base<P>(v);
    }

    template<Vector V>
    template<int Size>
    auto Conjugate<V>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return vec.template packet<Size>(index).conjugate();
    }

    template<Vector V>
    template<int Size>
    auto Conjugate<V>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return vec.template packet<Size>(index, count).conjugate();
    }

    template<Vector V>
    decltype(auto) Conjugate<V>::conjugate(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Vector V>
    __host__ __device__ consteval size_t Conjugate<V>::getSizeAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }
}

namespace Physica {
    template<Vector V>
    class Traits<Conjugate<V>> : public Traits<V> {
        static_assert(std::remove_cvref_t<V>::isComplex(), "[Error]: Unnecessary conjugate call on real vector");
    };
}

#ifdef PHYSICA_MKL
    #include "Conjugate_MKL.h"
#endif
