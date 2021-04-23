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
    class DenseMatrixType {
    public:
        enum {
            Column = 0b00,
            Row = 0b01,
            Vector = 0b00,
            Element = 0b10
        };
    public:
        template<class Matrix>
        constexpr static bool isColumnMatrix() { return !(Internal::Traits<Matrix>::MatrixType & Row); }

        template<class Matrix>
        constexpr static bool isRowMatrix() { return !isColumnMatrix<Matrix>(); }

        template<class Matrix>
        constexpr static bool isElementMatrix() { return Internal::Traits<Matrix>::MatrixType & Element; }

        template<class Matrix>
        constexpr static bool isVectorMatrix() { return !isElementMatrix<Matrix>(); }

        template<class Matrix1, class Matrix2>
        constexpr static bool isSameMajor() { return isColumnMatrix<Matrix1>() == isColumnMatrix<Matrix2>(); }
    private:
        DenseMatrixType();
    };
}