/*
 * Copyright 2022-2024 Weibo He.
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

namespace Physica::Core {
    template<class T>
    class RealMatrix : public RValueMatrix<RealMatrix<T>> {
        using This = RealMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& mat;
    public:
        RealMatrix(const T& mat_) : mat(mat_) {}
        RealMatrix(const This&) = delete;
        RealMatrix(This&&) noexcept = delete;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).real(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class ValueMatrix : public RValueMatrix<ValueMatrix<T>> {
        using This = ValueMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& mat;
    public:
        ValueMatrix(const T& mat_) : mat(mat_) {}
        ValueMatrix(const This&) = delete;
        ValueMatrix(This&&) noexcept = delete;
        ~ValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).value(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T, int GradOrder>
    class GradMatrix : public RValueMatrix<GradMatrix<T, GradOrder>> {
        using This = GradMatrix<T, GradOrder>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& mat;
    public:
        GradMatrix(const T& mat_) : mat(mat_) {}
        GradMatrix(const This&) = delete;
        GradMatrix(This&&) noexcept = delete;
        ~GradMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).template grad<GradOrder>(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };
}

namespace Physica {
    template <class T>
    class Traits<Core::RealMatrix<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template <class T>
    class Traits<Core::ValueMatrix<T>> {
    public:
        using ScalarType = T::ValueType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template <class T, int GradOrder>
    class Traits<Core::GradMatrix<T, GradOrder>> {
        static_assert(T::ScalarType::isDiffable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = T::ScalarType::template GradRtnTy<GradOrder>;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
