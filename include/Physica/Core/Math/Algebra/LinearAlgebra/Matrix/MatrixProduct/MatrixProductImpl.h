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
        constexpr static bool isAnyMajor = defaultMajor == MatrixOption::AnyMajor;
        using TargetType = LValueMatrix<OtherDerived>;

        if constexpr (isAnyMajor) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i) {
                for (size_t j = 0; j < target.getMaxMinor(); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<TargetType>(i, j);
                    const size_t c = MatrixOption::columnFromMajorMinor<TargetType>(i, j);
                    target.refFromMajorMinor(i, j) = calc(r, c);
                }
            }
        }
        else {
            for (size_t i = 0; i < (defaultMajor == MatrixOption::Column ? getColumn() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == MatrixOption::Column ?  getRow() : getColumn()); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<DefaultType>(i, j);
                    const size_t c = MatrixOption::columnFromMajorMinor<DefaultType>(i, j);
                    target(r, c) = calc(r, c);
                }
            }
        }

        constexpr bool isContinuous = std::is_base_of<ContinuousMatrix<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && isReverseDiff)
            target.getDerived().makeContinuous();
    }

    template<class T1, class T2>
    typename MatrixProduct<T1, T2>::ScalarType MatrixProduct<T1, T2>::calc(size_t row, size_t column) const {
        ScalarType result(0);
        for (size_t i = 0; i < mat1.getColumn(); ++i)
            result += ScalarType(mat1.calc(row, i)) * ScalarType(mat2.calc(i, column));
        return result;
    }

    template<class VectorType, class MatrixType>
    typename VectorMatrixProduct<VectorType, MatrixType>::ScalarType VectorMatrixProduct<VectorType, MatrixType>::calc(size_t row, size_t column) const {
        return vec.calc(row) * mat.calc(0, column);
    }

    template<class MatrixType, class VectorType>
    template<class OtherDerived, class Executor>
    inline void MatrixVectorProduct<MatrixType, VectorType>::assignTo(LValueVector<OtherDerived>& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);
        
        constexpr bool isContinuous = std::is_base_of<ContinuousVector<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && Base::isReverseDiff)
            target.getDerived().makeContinuous();
    }

    template<class MatrixType, class VectorType>
    inline typename MatrixVectorProduct<MatrixType, VectorType>::ScalarType MatrixVectorProduct<MatrixType, VectorType>::calc(size_t index) const {
        return mat.row(index) * vec;
    }
}
