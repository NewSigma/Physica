/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/ContinuousVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<class Derived> class RValueMatrix;
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class RMatrixBlock;

    template<Matrix T>
    class RMatrixBlock<T, 1, Dynamic> : public RValueVector<RMatrixBlock<T, 1, Dynamic>> {
        using This = RMatrixBlock<T, 1, Dynamic>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        T& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(T& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_) : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const {
            assert(index < colCount);
            return mat.calc(fromRow, fromCol + index);
        }
        [[nodiscard]] ValueType calc_value(size_t index) const {
            assert(index < colCount);
            return mat.calc_value(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<Matrix T>
    class RMatrixBlock<T, Dynamic, 1> : public RValueVector<RMatrixBlock<T, Dynamic, 1>> {
        using This = RMatrixBlock<T, Dynamic, 1>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        T& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        RMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_) : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const {
            assert(index < rowCount);
            return mat.calc(fromRow + index, fromCol);
        }
        [[nodiscard]] ValueType calc_value(size_t index) const {
            assert(index < rowCount);
            return mat.calc_value(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<Matrix T>
    class RMatrixBlock<T, Dynamic, Dynamic> : public RValueMatrix<RMatrixBlock<T, Dynamic, Dynamic>> {
        using This = RMatrixBlock<T, Dynamic, Dynamic>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        T& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
    };

    template<Matrix T>
    RMatrixBlock<T, Dynamic, Dynamic>::RMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix T>
    auto RMatrixBlock<T, Dynamic, Dynamic>::calc(size_t row, size_t col) const -> ScalarType {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }

    template<Matrix T>
    auto RMatrixBlock<T, Dynamic, Dynamic>::calc_value(size_t row, size_t col) const -> ValueType {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc_value(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    template<Matrix T, size_t Row, size_t Col>
    class Traits<Core::RMatrixBlock<T, Row, Col>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
