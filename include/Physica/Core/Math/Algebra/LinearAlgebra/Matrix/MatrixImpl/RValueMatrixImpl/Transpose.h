/*
 * Copyright 2021-2025 Weibo He.
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
    template<class T>
    struct remove_transpose {
        using Type = T;
    };

    template<Matrix T>
    struct remove_transpose<Transpose<T>> {
        using Type = T;
    };

    template<Vector T>
    struct remove_transpose<TransposeVector<T>> {
        using Type = T;
    };

    template<class T>
    struct is_transpose {
        constexpr static bool value = !std::is_same<T, typename remove_transpose<T>::Type>::value;
    };

    template<Matrix T>
    class Transpose<T> : public RValueMatrix<Transpose<T>> {
        using This = Transpose<T>;
        using Base = RValueMatrix<This>;

        const T& mat;
    public:        
        using typename Base::ScalarType;
        using typename Base::Tv;
    public:
        Transpose(const T& mat_) : mat(mat_) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(col, row); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return mat.calc_value(col, row); }
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
    };

    template<Vector T>
    class TransposeVector<T> : public RValueMatrix<TransposeVector<T>> {
        using This = TransposeVector<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::Tv;
    private:
        const T& vec;
    public:
        explicit TransposeVector(const T& vec_) : vec(vec_) {}
        TransposeVector(const This&) = default;
        TransposeVector(This&&) noexcept = default;
        ~TransposeVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] Tv calc_value([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc_value(col); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector T>
    void TransposeVector<T>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<Transpose<T>> {
    private:
        constexpr static int OtherMajor = MatrixOption::isColMatrix<T>() ? MatrixOption::Row : MatrixOption::Col;
        constexpr static int Major = MatrixOption::isAnyMajor<T>() ? MatrixOption::AnyMajor : OtherMajor;
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = Major;
        constexpr static size_t RowAtCompile = T::ColAtCompile;
        constexpr static size_t ColAtCompile = T::RowAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<Vector T>
    class Traits<TransposeVector<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = MatrixOption::Row;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColAtCompile = T::SizeAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
