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
    template<class OtherDerived>
    void MatrixProduct<MatrixType1, MatrixType2>::assignTo(LValueMatrix<OtherDerived>& target) const {
        constexpr static int defaultMajor = Internal::ProductOption<MatrixType1, MatrixType2>::Major;
        constexpr static bool isDynamic = defaultMajor == Dynamic;
        using TargetType = LValueMatrix<OtherDerived>;

        if constexpr (isDynamic) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i)
                for (size_t j = 0; j < target.getMaxMinor(); ++j)
                    target.getElementFromMajorMinor(i, j) = calc(TargetType::rowFromMajorMinor(i, j),
                                                                 TargetType::columnFromMajorMinor(i, j));
        }
        else {
            for (size_t i = 0; i < (defaultMajor == DenseMatrixOption::Column ? getColumn() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == DenseMatrixOption::Column ?  getRow() : getColumn()); ++j) {
                    const size_t r = DefaultType::rowFromMajorMinor(i, j);
                    const size_t c = DefaultType::columnFromMajorMinor(i, j);
                    target(r, c) = calc(r, c);
                }
            }
        }
    }

    template<class T1, class T2>
    typename MatrixProduct<T1, T2>::ScalarType MatrixProduct<T1, T2>::calc(size_t row, size_t column) const {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < mat1.getColumn(); ++i)
            result += ScalarType(mat1.calc(row, i) * mat2.calc(i, column));
        return result;
    }

    template<class VectorType, class MatrixType>
    template<class OtherDerived>
    void VectorMatrixProduct<VectorType, MatrixType>::assignTo(LValueMatrix<OtherDerived>& target) const {
        using TargetType = LValueMatrix<OtherDerived>;
        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                target.getElementFromMajorMinor(i, j) = calc(TargetType::rowFromMajorMinor(i, j),
                                                             TargetType::columnFromMajorMinor(i, j));
    }

    template<class VectorType, class MatrixType>
    typename VectorMatrixProduct<VectorType, MatrixType>::ScalarType VectorMatrixProduct<VectorType, MatrixType>::calc(size_t row, size_t column) const {
        return vec.calc(row) * mat.calc(0, column);
    }

    template<class MatrixType, class VectorType>
    template<class OtherDerived>
    void MatrixVectorProduct<MatrixType, VectorType>::assignTo(LValueVector<OtherDerived>& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);
    }

    template<class MatrixType, class VectorType>
    typename MatrixVectorProduct<MatrixType, VectorType>::ScalarType MatrixVectorProduct<MatrixType, VectorType>::calc(size_t index) const {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < vec.getLength(); ++i)
            result += ScalarType(mat.calc(index, i) * vec.calc(i));
        return result;
    }
}
