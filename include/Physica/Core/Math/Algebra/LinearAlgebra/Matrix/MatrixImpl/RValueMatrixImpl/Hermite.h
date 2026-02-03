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

    template<class T>
    struct is_hermite {
        constexpr static bool value = !std::is_same<T, typename remove_hermite<T>::Type>::value;
    };

    template<Matrix M>
    class Hermite<M> : public RValueMatrix<Hermite<M>> {
        using This = Hermite<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M& mat;
    public:
        explicit Hermite(const M& mat_);
        Hermite(const This&) = default;
        Hermite(This&&) noexcept = default;
        ~Hermite() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return mat.calc(col, row).conjugate(); }

        [[nodiscard]] const M& hermite() const noexcept { return mat; }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
    };

    template<Matrix M>
    Hermite<M>::Hermite(const M& mat_) : mat(mat_) {
        static_assert(M::isComplex, "[Error]: Do not call hermite on real matrix");
    }

    template<Vector V>
    class Hermite<V> : public RValueMatrix<Hermite<V>> {
        using This = Hermite<V>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const V& vec;
    public:
        explicit Hermite(const V& vec_);
        Hermite(const This&) = default;
        Hermite(This&&) noexcept = default;
        ~Hermite() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col).conjugate(); }

        [[nodiscard]] const V& hermite() const noexcept { return vec; }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    Hermite<V>::Hermite(const V& vec_) : vec(vec_) {
        static_assert(V::isComplex, "[Error]: Do not call hermite on real vector");
    }

    template<Vector V>
    void Hermite<V>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Hermite<M>> {
        constexpr static int OtherMajor = MatrixOption::isColMatrix<M>() ? MatrixOption::Row : MatrixOption::Col;
        constexpr static int Major = MatrixOption::isBothMajor<M>() ? MatrixOption::BothMajor : OtherMajor;
    public:
        using ScalarType = M::ScalarType;
        constexpr static int Option = Major;
        constexpr static size_t RowAtCompile = M::ColAtCompile;
        constexpr static size_t ColAtCompile = M::RowAtCompile;
        constexpr static size_t SizeAtCompile = M::SizeAtCompile;
    };

    template<Vector V>
    class Traits<Hermite<V>> {
    public:
        using ScalarType = V::ScalarType;
        constexpr static int Option = MatrixOption::Row;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColAtCompile = V::SizeAtCompile;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
    };
}
