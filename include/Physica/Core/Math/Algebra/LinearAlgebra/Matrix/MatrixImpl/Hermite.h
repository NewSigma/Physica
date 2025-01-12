/*
 * Copyright 2024 Weibo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    template<class MatrixType> class Hermite;
    template<class VectorType> class HermiteVector;

    template<class T>
    struct remove_hermite {
        using Type = T;
    };

    template<Matrix T>
    struct remove_hermite<Hermite<T>> {
        using Type = T;
    };

    template<Vector T>
    struct remove_hermite<HermiteVector<T>> {
        using Type = T;
    };

    template<class T>
    struct is_hermite {
        constexpr static bool value = std::is_same<T, typename remove_hermite<T>::Type>::value;
    };

    template<Matrix T>
    class Hermite<T> : public RValueMatrix<Hermite<T>> {
    public:
        using Base = RValueMatrix<Hermite<T>>;
        using typename Base::ScalarType;
    private:
        const T& matrix;
    public:
        Hermite(const T& matrix_) : matrix(matrix_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getRow(); }
    };

    template<Vector T>
    class HermiteVector<T> : public RValueMatrix<HermiteVector<T>> {
    public:
        using Base = RValueMatrix<HermiteVector<T>>;
        using typename Base::ScalarType;
    private:
        const T& vec;
    public:
        explicit HermiteVector(const T& vec_) : vec(vec_) {}
        /* Operations */
        template<Matrix M>
        void assignTo(M& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col).conjugate(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector T>
    template<Matrix M>
    void HermiteVector<T>::assignTo(M& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(M::rowFromMajorMinor(0, i), M::colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template <Matrix T>
    class Traits<Hermite<T>> {
    private:
        constexpr static int OtherMajor = MatrixOption::isColMatrix<T>() ? MatrixOption::Row : MatrixOption::Col;
        constexpr static int Major = MatrixOption::isAnyMajor<T>() ? MatrixOption::AnyMajor : OtherMajor;
        constexpr static int Storage = MatrixOption::getStorage<T>();
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = Major | Storage;
        constexpr static size_t RowAtCompile = T::ColAtCompile;
        constexpr static size_t ColAtCompile = T::RowAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template <Vector T>
    class Traits<HermiteVector<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColAtCompile = T::SizeAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
