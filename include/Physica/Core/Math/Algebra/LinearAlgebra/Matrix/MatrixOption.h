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

#include <cstdlib>
#include <utility>
#include <Physica/Macro.h>

namespace Physica::Core {
    /**
     * This enum decides how the data is stored in a matrix.
     * A dense matrix can be stored as elements or vectors in rows or columns.
     *
     * It is recommended that store elements to store a small matrix.
     */
    class MatrixOption {
    public:
        enum {
            Column = 0b0000,
            Row = 0b0001,
            AnyMajor = 0b0010,
            Vector = 0b0000,
            Element = 0b0100,
            AnyStorage = 0b1000
        };
    public:
        template<class Matrix>
        __host__ __device__ constexpr static bool isColumnMatrix() {
            return isAnyMajor<Matrix>() || !(Traits<Matrix>::Option & Row);
        }

        template<class Matrix>
        __host__ __device__ constexpr static bool isRowMatrix() {
            return isAnyMajor<Matrix>() || !isColumnMatrix<Matrix>();
        }

        template<class Matrix>
        __host__ __device__ constexpr static bool isAnyMajor() {
            return Traits<Matrix>::Option & AnyMajor;
        }

        template<class Matrix>
        __host__ __device__ constexpr static int getMajor() {
            return isAnyMajor<Matrix>() ? AnyMajor : (isColumnMatrix<Matrix>() ? Column : Row);
        }

        template<class Matrix>
        __host__ __device__ constexpr static bool isElementMatrix() {
            return isAnyStorage<Matrix>() || Traits<Matrix>::Option & Element;
        }

        template<class Matrix>
        __host__ __device__ constexpr static bool isVectorMatrix() {
            return isAnyStorage<Matrix>() || !isElementMatrix<Matrix>();
        }

        template<class Matrix>
        __host__ __device__ constexpr static bool isAnyStorage() {
            return Traits<Matrix>::Option & AnyStorage;
        }

        template<class Matrix>
        __host__ __device__ constexpr static int getStorage() {
            return isAnyStorage<Matrix>() ? AnyStorage : (isElementMatrix<Matrix>() ? Element : Vector);
        }

        template<class Matrix1, class Matrix2>
        __host__ __device__ constexpr static bool isSameMajor() { return isColumnMatrix<Matrix1>() == isColumnMatrix<Matrix2>(); }

        template<class Matrix1, class Matrix2>
        __host__ __device__ constexpr static bool isSameStorage() { return isElementMatrix<Matrix1>() == isElementMatrix<Matrix2>(); }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static size_t selectMajor(size_t row, size_t col) {
            return isColumnMatrix<Matrix>() ? col : row;
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static size_t selectMinor(size_t row, size_t col) {
            return isColumnMatrix<Matrix>() ? row : col;
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ static size_t getMaxMajor(const Matrix& mat) noexcept {
            return isColumnMatrix<Matrix>() ? mat.getColumn() : mat.getRow();
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ static size_t getMaxMinor(const Matrix& mat) noexcept {
            return isColumnMatrix<Matrix>() ? mat.getRow() : mat.getColumn();
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColumnMatrix<Matrix>() ? minor : major;
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static size_t columnFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColumnMatrix<Matrix>() ? major : minor;
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static bool isSymmMatrix() noexcept {
            using TransposeType = decltype(std::declval<Matrix>().transpose());
            using TransposeType1 = typename std::remove_cv<TransposeType>::type;
            using TransposeType2 = typename std::remove_reference<TransposeType1>::type;
            return std::is_base_of<TransposeType2, Matrix>::value;
        }

        template<class Matrix>
        [[nodiscard]] __host__ __device__ constexpr static bool isHermiteMatrix() noexcept {
            using HermiteType = decltype(std::declval<Matrix>().hermite());
            using HermiteType1 = typename std::remove_cv<HermiteType>::type;
            using HermiteType2 = typename std::remove_reference<HermiteType1>::type;
            return std::is_base_of<HermiteType2, Matrix>::value;
        }
    private:
        MatrixOption();
    };
}