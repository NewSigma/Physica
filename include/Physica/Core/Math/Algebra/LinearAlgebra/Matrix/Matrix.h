/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"

namespace Physica {
    template<class Derived> class RValueMatrix;

    namespace Internal {
        template<class T>
        concept MatrixObj = std::derived_from<T, RValueMatrix<T>>
                         || std::derived_from<T, device_obj<RValueMatrix<typename remove_device_obj<T>::type>>>;
    }

    template<class T>
    concept Matrix = Internal::MatrixObj<remove_codiff_t<std::remove_cvref_t<T>>>;

    class MatrixMajor {
    public:
        enum Option : char {
            Col = 0b01,
            Row = 0b10,
            BothMajor = 0b11, // Wildcard that indicates no performance difference between row access and col access

            RightMajor = Col,
            LeftMajor = Row,
            FMajor = Col,
            CMajor = Row,
        };
    public:
        MatrixMajor() = delete;
        /* Static members */
        template<class MatrixType>
        consteval static bool isColMatrix() noexcept {
            return (Traits<MatrixType>::Option & Col) != 0;
        }

        template<class MatrixType>
        consteval static bool isRowMatrix() noexcept {
            return (Traits<MatrixType>::Option & Row) != 0;
        }

        template<class MatrixType>
        consteval static bool isBothMajor() noexcept {
            return Traits<MatrixType>::Option == BothMajor;
        }

        template<class MatrixType>
        consteval static int getMajor() noexcept {
            return Traits<MatrixType>::Option;
        }

        template<class MatrixType1, class MatrixType2>
        consteval static bool isSameMajor() noexcept {
            return getMajor<MatrixType1>() == getMajor<MatrixType2>();
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t selectMajor(size_t row, size_t col) noexcept {
            return isColMatrix<MatrixType>() ? col : row;
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t selectMinor(size_t row, size_t col) noexcept {
            return isColMatrix<MatrixType>() ? row : col;
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t getMaxMajor(const MatrixType& mat) noexcept {
            return isColMatrix<MatrixType>() ? mat.getCol() : mat.getRow();
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t getMaxMinor(const MatrixType& mat) noexcept {
            return isColMatrix<MatrixType>() ? mat.getRow() : mat.getCol();
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<MatrixType>() ? minor : major;
        }

        template<class MatrixType>
        [[nodiscard]] constexpr static size_t colFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<MatrixType>() ? major : minor;
        }
    };
}