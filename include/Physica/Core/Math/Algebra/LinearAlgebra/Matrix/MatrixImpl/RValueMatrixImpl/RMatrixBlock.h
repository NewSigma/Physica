/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica {
    template<class Derived> class RValueMatrix;
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class RMatrixBlock;

    template<Matrix M>
    class RMatrixBlock<M, 1, Dynamic> : public RValueVector<RMatrixBlock<M, 1, Dynamic>> {
        using This = RMatrixBlock<M, 1, Dynamic>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount);
        RMatrixBlock(const This&) = default;
        RMatrixBlock(This&&) noexcept = default;
        ~RMatrixBlock() = default;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return colCount; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M>
    RMatrixBlock<M, 1, Dynamic>::RMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_)), fromRow(fromRow), fromCol(fromCol), colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, 1, Dynamic>::calc(size_t index) const -> T {
        assert(index < colCount);
        return mat.calc(fromRow, fromCol + index);
    }

    template<Matrix M>
    auto RMatrixBlock<M, 1, Dynamic>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return RMatrixBlock<M1, 1, Dynamic>(std::forward<M1>(v), self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M>
    class RMatrixBlock<M, Dynamic, 1> : public RValueVector<RMatrixBlock<M, Dynamic, 1>> {
        using This = RMatrixBlock<M, Dynamic, 1>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        RMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol);
        RMatrixBlock(const This&) = default;
        RMatrixBlock(This&&) noexcept = default;
        ~RMatrixBlock() = default;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return rowCount; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M>
    RMatrixBlock<M, Dynamic, 1>::RMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(std::forward<M>(mat_)), fromRow(fromRow), fromCol(fromCol), rowCount(rowCount) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, 1>::calc(size_t index) const -> T {
        assert(index < rowCount);
        return mat.calc(fromRow + index, fromCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, 1>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return RMatrixBlock<M1, Dynamic, 1>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol);
    }

    template<Matrix M>
    class RMatrixBlock<M, Dynamic, Dynamic> : public RValueMatrix<RMatrixBlock<M, Dynamic, Dynamic>> {
        using This = RMatrixBlock<M, Dynamic, Dynamic>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        RMatrixBlock(const This&) = default;
        RMatrixBlock(This&&) noexcept = default;
        ~RMatrixBlock() = default;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;
        [[nodiscard]] auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] auto topRows(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto bottomRows(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] auto leftCols(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto rightCols(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getCol() const noexcept { return colCount; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M>
    RMatrixBlock<M, Dynamic, Dynamic>::RMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        Base::checkBlock(mat, fromRow, rowCount, fromCol, colCount);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::calc(size_t row, size_t col) const -> T {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, 1, Dynamic>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::topRows(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::bottomRows(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, 0, self.getRow(), fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::leftCols(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, self.getRow(), 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::rightCols(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, 0, self.getRow(), from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        Base::checkBlock(self, 0, toRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::topLeftCorner(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        Base::checkBlock(self, 0, toRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    auto RMatrixBlock<M, Dynamic, Dynamic>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol, self.colCount);
    }

    template<Matrix M>
    __host__ __device__ consteval int RMatrixBlock<M, Dynamic, Dynamic>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<RMatrixBlock<M, Row, Col>> {
    public:
        using ScalarType = std::remove_cvref_t<M>::ScalarType;
    };
}
