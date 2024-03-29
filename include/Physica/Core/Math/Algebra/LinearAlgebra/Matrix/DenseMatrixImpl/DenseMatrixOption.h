/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits;
    }
    /**
     * This enum decides how the data is stored in a matrix.
     * A dense matrix can be stored as elements or vectors in rows or columns.
     *
     * It is recommended that store elements to store a small matrix.
     */
    class DenseMatrixOption {
    public:
        enum {
            Column = 0b000,
            Row = 0b001,
            AnyMajor = 0b010,
            Vector = 0b000,
            Element = 0b100
        };
    public:
        template<class Matrix>
        constexpr static bool isColumnMatrix() {
            return isAnyMajor<Matrix>() || !(Internal::Traits<Matrix>::MatrixOption & Row);
        }

        template<class Matrix>
        constexpr static bool isRowMatrix() {
            return isAnyMajor<Matrix>() || !isColumnMatrix<Matrix>();
        }

        template<class Matrix>
        constexpr static bool isAnyMajor() {
            return Internal::Traits<Matrix>::MatrixOption & AnyMajor;
        }

        template<class Matrix>
        constexpr static int getMajor() {
            return isAnyMajor<Matrix>() ? AnyMajor : (isColumnMatrix<Matrix>() ? Column : Row);
        }

        template<class Matrix>
        constexpr static bool isElementMatrix() { return Internal::Traits<Matrix>::MatrixOption & Element; }

        template<class Matrix>
        constexpr static bool isVectorMatrix() { return !isElementMatrix<Matrix>(); }

        template<class Matrix>
        constexpr static int getStorage() { return isElementMatrix<Matrix>() ? Element : Vector; }

        template<class Matrix1, class Matrix2>
        constexpr static bool isSameMajor() { return isColumnMatrix<Matrix1>() == isColumnMatrix<Matrix2>(); }
    private:
        DenseMatrixOption();
    };
}