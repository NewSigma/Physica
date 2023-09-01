/*
 * Copyright 2023 WeiBo He.
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
    template<class Derived>
    inline device_obj<ContinuousMatrix<Derived>>&
    device_obj<ContinuousMatrix<Derived>>::operator=(const device_obj<ContinuousMatrix<Derived>>& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<class Derived>
    inline device_obj<ContinuousMatrix<Derived>>&
    device_obj<ContinuousMatrix<Derived>>::operator=(device_obj<ContinuousMatrix<Derived>>&& obj) noexcept {
        Base::operator=(std::forward<Base>(obj));
        return *this;
    }

    template<class Derived>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, 1, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::row(size_t r) {
        const bool useSpecialization = device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getDerived(), r, 1, 0};
        else
            return {Base::getDerived(), r, 0, Base::getColumn()};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, 1, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::row(size_t r) const {
        const bool useSpecialization = device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getConstCastDerived(), r, 1, 0};
        else
            return {Base::getConstCastDerived(), r, 0, Base::getColumn()};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, 1>>
    device_obj<ContinuousMatrix<Derived>>::col(size_t c) {
        return {Base::getDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, 1>>
    device_obj<ContinuousMatrix<Derived>>::col(size_t c) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::bottomRows(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, device_obj<ContinuousMatrix<Derived>>::ColumnAtCompile>>
    device_obj<ContinuousMatrix<Derived>>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::leftCols(size_t to) {
        return {Base::getDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::rightCols(size_t from) {
        return {Base::getDerived(), 0, Base::getRow(), from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, device_obj<ContinuousMatrix<Derived>>::RowAtCompile, Column>>
    device_obj<ContinuousMatrix<Derived>>::rightCols(size_t from) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) {
        return {Base::getDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) const {
        return {Base::getConstCastDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t to) {
        return {Base::getDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) {
        return {Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return {Base::getDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>>
    device_obj<ContinuousMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount};
    }
}
