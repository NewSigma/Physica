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
    template<class MatrixType1, class MatrixType2>
    class MatrixProduct {
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename MatrixType1::ScalarType,
                                                                       typename MatrixType2::ScalarType>::Type;

        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixProduct(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t column) const;
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat2.getColumn(); }
    };

    template<class Derived, class OtherDerived>
    inline MatrixProduct<Derived, OtherDerived>
    operator*(const RValueMatrix<Derived>& mat1, const RValueMatrix<OtherDerived>& mat2) {
        assert(mat1.getColumn() == mat2.getRow());
        return MatrixProduct(mat1, mat2);
    }

    template<class T1, class T2>
    typename MatrixProduct<T1, T2>::ScalarType MatrixProduct<T1, T2>::operator()(size_t row, size_t column) const {
        ScalarType result = 0;
        for (size_t i = 0; i < mat1.getColumn(); ++i)
            result += ScalarType(mat1(row, i) * mat2(i, column));
        return result;
    }
}
