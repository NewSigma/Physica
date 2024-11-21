/*
 * Copyright 2021-2024 Weibo He.
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
    template<Matrix T>
    class Transpose<T> : public RValueMatrix<Transpose<T>> {
        using Base = RValueMatrix<Transpose<T>>;

        const T& matrix;
    public:        
        using typename Base::ScalarType;
    public:
        Transpose(const T& matrix_) : matrix(matrix_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getRow(); }
    };

    template<Vector T>
    class TransposeVector<T> : public RValueMatrix<TransposeVector<T>> {
        using This = TransposeVector<T>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const T& vec;
    public:
        explicit TransposeVector(const T& vec_) : vec(vec_) {}
        /* Operations */
        template<Matrix M>
        void assignTo(LValueMatrix<M>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector T>
    template<Matrix M>
    void TransposeVector<T>::assignTo(LValueMatrix<M>& target) const {
        using TargetType = LValueMatrix<M>;
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(TargetType::rowFromMajorMinor(0, i), TargetType::colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<Transpose<T>> {
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

    template<Vector T>
    class Traits<TransposeVector<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColAtCompile = T::SizeAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
