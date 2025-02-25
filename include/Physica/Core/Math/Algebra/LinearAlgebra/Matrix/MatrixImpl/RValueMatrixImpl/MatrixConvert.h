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

namespace Physica {
    template<class T>
    class RealMatrix : public RValueMatrix<RealMatrix<T>> {
        using This = RealMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& mat;
    public:
        RealMatrix(const T& mat_) : mat(mat_) {}
        RealMatrix(const This&) = default;
        RealMatrix(This&&) noexcept = default;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).real(); }
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class ImagMatrix : public RValueMatrix<ImagMatrix<T>> {
        using This = ImagMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& mat;
    public:
        ImagMatrix(const T& mat_) : mat(mat_) {}
        ImagMatrix(const This&) = default;
        ImagMatrix(This&&) noexcept = default;
        ~ImagMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).imag(); }
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class SquaredNormMatrix : public RValueMatrix<SquaredNormMatrix<T>> {
        using This = SquaredNormMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& mat;
    public:
        SquaredNormMatrix(const T& mat_) : mat(mat_) {}
        SquaredNormMatrix(const This&) = default;
        SquaredNormMatrix(This&&) noexcept = default;
        ~SquaredNormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).squaredNorm(); }
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().squaredNorm(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class NormMatrix : public RValueMatrix<NormMatrix<T>> {
        using This = NormMatrix<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& mat;
    public:
        NormMatrix(const T& mat_) : mat(mat_) {}
        NormMatrix(const This&) = default;
        NormMatrix(This&&) noexcept = default;
        ~NormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).norm(); }
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().norm(); }
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
        ValueMatrix(const This&) = default;
        ValueMatrix(This&&) noexcept = default;
        ~ValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc_value(row, col); }
        [[nodiscard]] ScalarType calc_value(size_t row, size_t col) const { return calc(row, col); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class T, int GradOrder>
    class GradMatrix : public RValueMatrix<GradMatrix<T, GradOrder>> {
        using This = GradMatrix<T, GradOrder>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& mat;
    public:
        GradMatrix(const T& mat_) : mat(mat_) {}
        GradMatrix(const This&) = default;
        GradMatrix(This&&) noexcept = default;
        ~GradMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).template grad<GradOrder>(); }
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<RealMatrix<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<class T>
    class Traits<ImagMatrix<T>> : public Traits<RealMatrix<T>> {};

    template<class T>
    class Traits<SquaredNormMatrix<T>> : public Traits<RealMatrix<T>> {};

    template<class T>
    class Traits<NormMatrix<T>> : public Traits<RealMatrix<T>> {};

    template<class T>
    class Traits<ValueMatrix<T>> {
    public:
        using ScalarType = T::ValueType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<class T, int GradOrder>
    class Traits<GradMatrix<T, GradOrder>> {
        static_assert(T::ScalarType::isDiffable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = Internal::GradTypeHelper<typename T::ScalarType, GradOrder>::Type;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
