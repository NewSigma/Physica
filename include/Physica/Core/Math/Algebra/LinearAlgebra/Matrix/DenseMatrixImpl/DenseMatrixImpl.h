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

#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica::Core {
    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(size_t row, size_t column) : Storage(row, column), Dim(row, column) {}

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(size_t row, size_t column, T value) : Storage(row, column, std::move(value)), Dim(row, column) {}

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(std::initializer_list<InitializerType> list) : Storage(std::move(list)) {
        if constexpr (MatrixOption::isVectorMatrix<This>()) {
            const size_t length = Storage::asArray().getLength();
            const size_t size = Storage::getSize();
            if constexpr (MatrixOption::isColumnMatrix<This>())
                Dim::resize(size / length, length);
            else
                Dim::resize(length, size / length);
        }
        else {
            Dim::resize(Row, Column);
        }
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class OtherMatrix>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(const RValueMatrix<OtherMatrix>& mat)
            : DenseMatrix(mat.getRow(), mat.getColumn()) {
        mat.getDerived().assignTo(*this);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class VectorType>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(const RValueVector<VectorType>& vec)
            : DenseMatrix(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.getDerived().assignTo(col);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class MatrixIn>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(LUDecomposition<MatrixIn> lu)
            : DenseMatrix(lu.getRow(), lu.getRow()) {
        const size_t rank = lu.getRow();
        (*this) = lu.getMatrix();
        for (size_t i = 0; i < rank; ++i)
            lu.decompositionColumn((*this), i);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(const DenseMatrix& m)
            : Base(), Storage(static_cast<const Storage&>(m)), Dim(m) {}

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::DenseMatrix(DenseMatrix&& m) noexcept : Base(), Storage(std::move(m)), Dim(m) {}

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>&
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::operator=(DenseMatrix m) noexcept {
        swap(m);
        return *this;
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::resize(size_t row, size_t column) {
        Storage::resize(row, column);
        Dim::resize(row, column);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::swap(DenseMatrix& m) noexcept {
        Storage::swap(m);
        Dim::swap(m);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class VectorType>
    std::pair<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>, DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::meshgrid(
            const LValueVector<VectorType>& vecX,
            const LValueVector<VectorType>& vecY) {
        using MatrixType = DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>;
        const size_t row = vecY.getLength();
        const size_t col = vecX.getLength();
        MatrixType x(row, col);
        MatrixType y(row, col);
        for (size_t i = 0; i < x.getMaxMajor(); ++i) {
            for (size_t j = 0; j < x.getMaxMinor(); ++j) {
                x.refFromMajorMinor(i, j) = vecX[MatrixOption::columnFromMajorMinor<MatrixType>(i, j)];
                y.refFromMajorMinor(i, j) = vecY[MatrixOption::rowFromMajorMinor<MatrixType>(i, j)];
            }
        }
        return std::make_pair(std::move(x), std::move(y));
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator> DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::unitMatrix(size_t order) {
        DenseMatrix result(order, order);
        result.toUnitMatrix();
        return result;
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class RandomGenerator>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator> DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        DenseMatrix result(row, column);
        for (size_t i = 0; i < result.getMaxMajor(); ++i)
            for (size_t j = 0; j < result.getMaxMinor(); ++j)
                result.refFromMajorMinor(i, j) = T::random_uniform(gen);
        return result;
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>& mat) {
        const size_t column = mat.getColumn();
        const size_t row = mat.getRow();
        size_t width = 0;
        /* Get max width */ {
            for (size_t c = 0; c < column; ++c) {
                for (size_t r = 0; r < row; ++ r) {
                    std::stringstream stream{};
                    stream.copyfmt(os);
                    stream << mat(r, c).getReal();
                    width = std::max(width, stream.str().length());
                }
            }
        }
        /* Output */ {
            for (size_t r = 0; r < row; ++r) {
                for (size_t c = 0; c < column; ++c) {
                    os.width(width);
                    os << mat(r, c) << ' ';
                }
                os << '\n';
            }
        }
        return os;
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    std::istream& operator>>(std::istream& is, DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>& mat) {
        const size_t column = mat.getColumn();
        const size_t row = mat.getRow();
        for (size_t r = 0; r < row; ++r)
            for (size_t c = 0; c < column; ++c)
                is >> mat(r, c);
        if (!is)
            throw BadFileFormatException();
        return is;
    }
}
