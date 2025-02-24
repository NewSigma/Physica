/*
 * Copyright 2023 Weibo He.
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

#include "../ContinuousMatrix.cuh"

namespace Physica {
    template<class Derived>
    inline auto device_obj<ContinuousMatrix<Derived>>::operator=(const This& obj) -> This& {
        Base::operator=(obj);
        return *this;
    }

    template<class Derived>
    inline auto device_obj<ContinuousMatrix<Derived>>::operator=(This&& obj) noexcept -> This& {
        Base::operator=(std::forward<Base>(obj));
        return *this;
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::row(size_t r) noexcept {
        const bool useSpecialization = device_obj<ContinuousMatrix<Derived>>::ColAtCompile == 1;
        if constexpr (useSpecialization)
            return RowVector(Base::getDerived(), r, 1, 0);
        else
            return RowVector(Base::getDerived(), r, 0, Base::getCol());
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::row(size_t r) const noexcept {
        const bool useSpecialization = device_obj<ContinuousMatrix<Derived>>::ColAtCompile == 1;
        if constexpr (useSpecialization)
            return RowVector(Base::getConstCastDerived(), r, 1, 0);
        else
            return RowVector(Base::getConstCastDerived(), r, 0, Base::getCol());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) noexcept {
        return RowBlock<Row>(Base::getDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::topRows(size_t to) noexcept {
        return RowBlock<Row>(Base::getDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::topRows(size_t to) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::bottomRows(size_t from) noexcept {
        return RowBlock<Row>(Base::getDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::bottomRows(size_t from) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::cols(size_t fromCol, size_t colCount) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::cols(size_t fromCol, size_t colCount) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::leftCols(size_t to) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::leftCols(size_t to) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::rightCols(size_t from) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::rightCols(size_t from) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t to) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::topLeftCorner(size_t to) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t from) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::bottomRightCorner(size_t from) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<ContinuousMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ inline const auto device_obj<ContinuousMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    auto device_obj<ContinuousMatrix<Derived>>::flatten() {
        return device_obj<FlattenC<Derived>>(Base::getDerived());
    }

    template<class Derived>
    const auto device_obj<ContinuousMatrix<Derived>>::flatten() const {
        return device_obj<FlattenC<Derived>>(const_cast<This&>(*this));
    }
}
