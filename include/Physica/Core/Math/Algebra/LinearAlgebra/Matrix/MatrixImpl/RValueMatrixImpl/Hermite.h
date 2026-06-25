/*
 * Copyright 2024-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<class> class Hermite;

    template<class T>
    struct remove_hermite {
        using Type = T;
    };

    template<class T>
    struct remove_hermite<Hermite<T>> {
        using Type = T;
    };

    template<Matrix M>
    class Hermite<M> : public RValueMatrix<Hermite<M>> {
        using This = Hermite<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> m;
    public:
        explicit Hermite(M&& m);
        Hermite(const This&) = default;
        Hermite(This&&) noexcept = default;
        ~Hermite() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return m.calc(col, row).conjugate(); }

        [[nodiscard]] auto&& hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return m.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return m.getRow(); }
        [[nodiscard]] size_t getOrder() const noexcept { return m.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M>
    Hermite<M>::Hermite(M&& m) : m(std::forward<M>(m)) {
        static_assert(m.isComplex(), "[Error]: Do not call hermite on real matrix");
    }

    template<Matrix M>
    auto&& Hermite<M>::hermite(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M>
    __host__ __device__ consteval bool Hermite<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Hermite<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Hermite<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval int Hermite<M>::getMajor() noexcept {
        return Transpose<M>::getMajor();
    }

    template<Vector V>
    class Hermite<V> : public RValueMatrix<Hermite<V>> {
        using This = Hermite<V>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<V> v;
    public:
        explicit Hermite(V&& v);
        Hermite(const This&) = default;
        Hermite(This&&) noexcept = default;
        ~Hermite() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return v.calc(col).conjugate(); }

        [[nodiscard]] auto&& hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getCol() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Vector V>
    Hermite<V>::Hermite(V&& v) : v(std::forward<V>(v)) {
        static_assert(v.isComplex(), "[Error]: Do not call hermite on real vector");
    }

    template<Vector V>
    void Hermite<V>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < v.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }

    template<Vector V>
    auto&& Hermite<V>::hermite(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v);
    }

    template<Vector V>
    __host__ __device__ consteval size_t Hermite<V>::getRowAtCompile() noexcept {
        return 1;
    }

    template<Vector V>
    __host__ __device__ consteval size_t Hermite<V>::getColAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }

    template<Vector V>
    __host__ __device__ consteval int Hermite<V>::getMajor() noexcept {
        return MatrixMajor::Row;
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Hermite<M>> : public Traits<Transpose<M>> {};

    template<Vector V>
    class Traits<Hermite<V>> : public Traits<Transpose<V>> {};
}
