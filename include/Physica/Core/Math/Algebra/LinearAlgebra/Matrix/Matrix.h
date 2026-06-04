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
            FMajor = Col, // Fortran convention
            CMajor = Row, // C conversion
        };
    public:
        MatrixMajor() = delete;
        /* Static members */
        template<class M>
        consteval static bool isColMatrix() noexcept {
            return (getMajor<M>() & Col) != 0;
        }

        template<class M>
        consteval static bool isRowMatrix() noexcept {
            return (getMajor<M>() & Row) != 0;
        }

        template<class M>
        consteval static bool isBothMajor() noexcept {
            return getMajor<M>() == BothMajor;
        }

        template<class M>
        consteval static int getMajor() noexcept {
            return std::remove_cvref_t<M>::getMajor();
        }

        template<class M1, class M2>
        consteval static bool isSameMajor() noexcept {
            return getMajor<M1>() == getMajor<M2>();
        }

        template<class M>
        [[nodiscard]] constexpr static size_t selectMajor(size_t row, size_t col) noexcept {
            return isColMatrix<M>() ? col : row;
        }

        template<class M>
        [[nodiscard]] constexpr static size_t selectMinor(size_t row, size_t col) noexcept {
            return isColMatrix<M>() ? row : col;
        }

        [[nodiscard]] constexpr static size_t getMaxMajor(const Matrix auto& mat) noexcept {
            return isColMatrix<decltype(mat)>() ? mat.getCol() : mat.getRow();
        }

        [[nodiscard]] constexpr static size_t getMaxMinor(const Matrix auto& mat) noexcept {
            return isColMatrix<decltype(mat)>() ? mat.getRow() : mat.getCol();
        }

        template<class M>
        [[nodiscard]] constexpr static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<M>() ? minor : major;
        }

        template<class M>
        [[nodiscard]] constexpr static size_t colFromMajorMinor(size_t major, size_t minor) noexcept {
            return isColMatrix<M>() ? major : minor;
        }
    };
}