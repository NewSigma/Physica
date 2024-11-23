/*
 * Copyright 2021 Weibo He.
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

#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica::Core {
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::DenseMatrix(size_t row, size_t col)
            : Storage(row, col) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::DenseMatrix(size_t row, size_t col, T value)
            : Storage(row, col, std::move(value)) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::DenseMatrix(std::initializer_list<InitializerType> list)
            : Storage(std::move(list)) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<Matrix M>
    DenseMatrix<T, Option, Row, Col, Allocator>::DenseMatrix(const M& mat)
            : DenseMatrix(mat.getRow(), mat.getCol()) {
        mat.assignTo(*this);
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<Vector V>
    DenseMatrix<T, Option, Row, Col, Allocator>::DenseMatrix(const V& vec)
            : DenseMatrix(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assignTo(col);
    }
    /**
     * \returns the origin column index of the main element
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:35
     */
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    size_t DenseMatrix<T, Option, Row, Col, Allocator>::completePivoting(size_t col) {
        const auto rank = getRow();
        assert(col < rank);
        size_t main_row_index = 0, main_col_index = 0;
        const T zero = T(0);
        const T* main = &zero;
        for (size_t i = col; i < rank; ++i) {
            for (size_t j = col; j < rank; ++j) {
                const auto* temp = Base::data_ptr(i, j);
                bool larger = absCompare(*main, *temp);
                main = larger ? main : temp;
                main_row_index = larger ? main_row_index : j;
                main_col_index = larger ? main_col_index : i;
            }
        }

        if (col != main_row_index)
            Storage::swap_row(col, main_row_index);
        return main_col_index;
    }
    /**
     * \returns the origin col index of the main element
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:35
     */
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    size_t DenseMatrix<T, Option, Row, Col, Allocator>::partialPivoting(size_t col) {
        const auto rank = getRow();
        assert(col < rank);
        size_t main_col_index = col;
        const T* main = Base::data_ptr(col, col);
        for (size_t j = col + 1; j < rank; ++j) {
            const auto* temp = Base::data_ptr(j, col);
            bool larger = absCompare(*main, *temp);
            main = larger ? main : temp;
            main_col_index = larger ? main_col_index : j;
        }

        if (col != main_col_index)
            Storage::swap_row(col, main_col_index);
        return main_col_index;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::copy() const {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            TracerType::getInstance().reserve(Storage::getSize());
            This result(Base::getRow(), Base::getCol());
            for (size_t major = 0; major < Base::getMaxMajor(); ++major)
                for (size_t minor = 0; minor < Base::getMaxMinor(); ++minor)
                    result.refFromMajorMinor(major, minor) = Base::refFromMajorMinor(major, minor).copy();
            return result;
        }
        else
            return *this;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator> DenseMatrix<T, Option, Row, Col, Allocator>::unitMatrix(size_t order) {
        DenseMatrix result(order, order);
        result.toUnitMatrix();
        return result;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<RandomGenerator R>
    inline DenseMatrix<T, Option, Row, Col, Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::random_uniform(
            size_t row, size_t col) {
        DenseMatrix result(row, col);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<RandomGenerator R>
    inline DenseMatrix<T, Option, Row, Col, Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::random_normal(
            size_t row, size_t col) {
        DenseMatrix result(row, col);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<class Distribution, RandomGenerator R>
    inline DenseMatrix<T, Option, Row, Col, Allocator>
    DenseMatrix<T, Option, Row, Col, Allocator>::random_any(
            size_t row, size_t col, Distribution& dist) {
        DenseMatrix result(row, col);
        result.template random_any<R>(dist);
        return result;
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<Vector V>
    std::pair<DenseMatrix<T, Option, Row, Col, Allocator>, DenseMatrix<T, Option, Row, Col, Allocator>>
    DenseMatrix<T, Option, Row, Col, Allocator>::meshgrid(const V& vecX, const V& vecY) {
        using MatrixType = DenseMatrix<T, Option, Row, Col, Allocator>;
        const size_t row = vecY.getLength();
        const size_t col = vecX.getLength();
        MatrixType x(row, col);
        MatrixType y(row, col);
        for (size_t i = 0; i < x.getMaxMajor(); ++i) {
            for (size_t j = 0; j < x.getMaxMinor(); ++j) {
                x.refFromMajorMinor(i, j) = vecX.calc(MatrixOption::colFromMajorMinor<MatrixType>(i, j));
                y.refFromMajorMinor(i, j) = vecY.calc(MatrixOption::rowFromMajorMinor<MatrixType>(i, j));
            }
        }
        return std::make_pair(std::move(x), std::move(y));
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    std::istream& operator>>(std::istream& is, DenseMatrix<T, Option, Row, Col, Allocator>& mat) {
        const size_t col = mat.getCol();
        const size_t row = mat.getRow();
        for (size_t r = 0; r < row; ++r)
            for (size_t c = 0; c < col; ++c)
                is >> mat(r, c);
        if (!is)
            throw BadFileFormatException("[Error]: bad matrix format");
        return is;
    }
}
