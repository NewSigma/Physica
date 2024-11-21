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

#include <concepts>
#include <cstdlib>
#include <utility>
#include "Physica/Macro.h"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    template<class Derived> class RValueMatrix;
    template<class Derived> class LValueMatrix;

    template<class T>
    concept Matrix = std::derived_from<T, RValueMatrix<T>> || std::derived_from<T, device_obj<RValueMatrix<T>>>;

    template<class T>
    concept LMatrix = std::derived_from<T, LValueMatrix<T>> || std::derived_from<T, device_obj<LValueMatrix<T>>>;
    /**
     * This enum decides how the data is stored in a matrix.
     * A dense matrix can be stored as elements or vectors in rows or cols.
     *
     * It is recommended that store elements to store a small matrix.
     */
    class MatrixOption {
    public:
        enum {
            Col = 0b0000,
            Row = 0b0001,
            AnyMajor = 0b0010,
            Vector = 0b0000,
            Element = 0b0100,
            AnyStorage = 0b1000
        };
    public:
        template<class MatrixType>
        __host__ __device__ constexpr static bool isColMatrix() {
            return isAnyMajor<MatrixType>() || !(Traits<MatrixType>::Option & Row);
        }

        template<class MatrixType>
        __host__ __device__ constexpr static bool isRowMatrix() {
            return isAnyMajor<MatrixType>() || !isColMatrix<MatrixType>();
        }

        template<class MatrixType>
        __host__ __device__ constexpr static bool isAnyMajor() {
            return Traits<MatrixType>::Option & AnyMajor;
        }

        template<class MatrixType>
        __host__ __device__ constexpr static int getMajor() {
            return isAnyMajor<MatrixType>() ? AnyMajor : (isColMatrix<MatrixType>() ? Col : Row);
        }

        template<class MatrixType>
        __host__ __device__ constexpr static bool isElementMatrix() {
            return isAnyStorage<MatrixType>() || Traits<MatrixType>::Option & Element;
        }

        template<class MatrixType>
        __host__ __device__ constexpr static bool isVectorMatrix() {
            return isAnyStorage<MatrixType>() || !isElementMatrix<MatrixType>();
        }

        template<class MatrixType>
        __host__ __device__ constexpr static bool isAnyStorage() {
            return Traits<MatrixType>::Option & AnyStorage;
        }

        template<class MatrixType>
        __host__ __device__ constexpr static int getStorage() {
            return isAnyStorage<MatrixType>() ? AnyStorage : (isElementMatrix<MatrixType>() ? Element : Vector);
        }

        template<class MatrixType1, class MatrixType2>
        __host__ __device__ constexpr static bool isSameMajor() { return getMajor<MatrixType1>() == getMajor<MatrixType2>(); }

        template<class MatrixType1, class MatrixType2>
        __host__ __device__ constexpr static bool isSameStorage() { return getStorage<MatrixType1>() == getStorage<MatrixType2>(); }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static size_t selectMajor(size_t row, size_t col) {
            return isColMatrix<MatrixType>() ? col : row;
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static size_t selectMinor(size_t row, size_t col) {
            return isColMatrix<MatrixType>() ? row : col;
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ static size_t getMaxMajor(const MatrixType& mat) noexcept {
            return isColMatrix<MatrixType>() ? mat.getCol() : mat.getRow();
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ static size_t getMaxMinor(const MatrixType& mat) noexcept {
            return isColMatrix<MatrixType>() ? mat.getRow() : mat.getCol();
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<MatrixType>() ? minor : major;
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static size_t colFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<MatrixType>() ? major : minor;
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static bool isSymmMatrix() noexcept {
            using TransposeType = decltype(std::declval<MatrixType>().transpose());
            using TransposeType1 = std::remove_cv<TransposeType>::type;
            using TransposeType2 = std::remove_reference<TransposeType1>::type;
            return std::is_base_of<TransposeType2, MatrixType>::value;
        }

        template<class MatrixType>
        [[nodiscard]] __host__ __device__ constexpr static bool isHermiteMatrix() noexcept {
            using HermiteType = decltype(std::declval<MatrixType>().hermite());
            using HermiteType1 = std::remove_cv<HermiteType>::type;
            using HermiteType2 = std::remove_reference<HermiteType1>::type;
            return std::is_base_of<HermiteType2, MatrixType>::value;
        }
    private:
        MatrixOption();
    };
}