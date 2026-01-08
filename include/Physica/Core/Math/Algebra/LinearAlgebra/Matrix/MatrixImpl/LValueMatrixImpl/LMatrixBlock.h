/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"
#include "../LValueMatrix.h"

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class LMatrixBlock;

    template<Matrix M, size_t Col>
    class LMatrixBlock<M, 1, Col> : public LValueVector<LMatrixBlock<M, 1, Col>> {
        using This = LMatrixBlock<M, 1, Col>;
        using Base = LValueVector<This>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        M& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(M& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_) : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Matrix M, size_t Col>
    auto LMatrixBlock<M, 1, Col>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.colCount);
        return self.mat.data_ptr(self.fromRow, self.fromCol + index);
    }

    template<Matrix M, size_t Row>
    class LMatrixBlock<M, Row, 1> : public LValueVector<LMatrixBlock<M, Row, 1>> {
        using This = LMatrixBlock<M, Row, 1>;
        using Base = LValueVector<This>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        M& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(M& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Matrix M, size_t Row>
    auto LMatrixBlock<M, Row, 1>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.rowCount);
        return self.mat.data_ptr(self.fromRow + index, self.fromCol);
    }

    template<Matrix M>
    class LMatrixBlock<M, 1, 1> : public LValueVector<LMatrixBlock<M, 1, 1>> {
        using This = LMatrixBlock<M, 1, 1>;
        using Base = LValueVector<This>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        M& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(M& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] auto data_ptr(this auto&& self, [[maybe_unused]] size_t index) noexcept;
    };

    template<Matrix M>
    auto LMatrixBlock<M, 1, 1>::data_ptr(this auto&& self, [[maybe_unused]] size_t index) noexcept {
        assert(index == 0 && "[Error]: Index overflow");
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    class LMatrixBlock<M, Dynamic, Dynamic> : public LValueMatrix<LMatrixBlock<M, Dynamic, Dynamic>> {
        using This = LMatrixBlock<M, Dynamic, Dynamic>;
        using Base = LValueMatrix<This>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        M& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(M& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getCol() const noexcept { return colCount; }
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
    };

    template<Matrix M>
    LMatrixBlock<M, Dynamic, Dynamic>::LMatrixBlock(M& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.rowCount);
        assert(col < self.colCount);
        return self.mat.data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<LMatrixBlock<M, Row, Col>> {
    public:
        using ScalarType = M::ScalarType;
        constexpr static int Option = M::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
