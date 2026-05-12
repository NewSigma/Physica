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
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount);
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Col; }
    };

    template<Matrix M, size_t Col>
    LMatrixBlock<M, 1, Col>::LMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_)), fromRow(fromRow), fromCol(fromCol), colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M, size_t Col>
    auto LMatrixBlock<M, 1, Col>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.colCount);
        return self.mat.data_ptr(self.fromRow, self.fromCol + index);
    }

    template<Matrix M, size_t Col>
    auto LMatrixBlock<M, 1, Col>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        if constexpr (v.isLValueMatrix())
            return LMatrixBlock<M1, 1, Col>(std::forward<M1>(v), self.fromRow, self.fromCol, self.getLength());
        else
            return RMatrixBlock<M1, 1, Dynamic>(std::forward<M1>(v), self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M, size_t Row>
    class LMatrixBlock<M, Row, 1> : public LValueVector<LMatrixBlock<M, Row, 1>> {
        using This = LMatrixBlock<M, Row, 1>;
        using Base = LValueVector<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol);
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Row; }
    };

    template<Matrix M, size_t Row>
    LMatrixBlock<M, Row, 1>::LMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(std::forward<M>(mat_)), fromRow(fromRow), fromCol(fromCol), rowCount(rowCount) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M, size_t Row>
    auto LMatrixBlock<M, Row, 1>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.rowCount);
        return self.mat.data_ptr(self.fromRow + index, self.fromCol);
    }

    template<Matrix M, size_t Row>
    auto LMatrixBlock<M, Row, 1>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        if constexpr (v.isLValueMatrix())
            return LMatrixBlock<M1, Row, 1>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol);
        else
            return RMatrixBlock<M1, Dynamic, 1>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol);
    }

    template<Matrix M>
    class LMatrixBlock<M, 1, 1> : public LValueVector<LMatrixBlock<M, 1, 1>> {
        using This = LMatrixBlock<M, 1, 1>;
        using Base = LValueVector<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t fromCol;
    public:
        LMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol);
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == 1); }
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] auto data_ptr(this auto&& self, [[maybe_unused]] size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return 1; }
    };

    template<Matrix M>
    LMatrixBlock<M, 1, 1>::LMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol)
            : mat(std::forward<M>(mat_)), fromRow(fromRow), fromCol(fromCol) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, 1, 1>::data_ptr(this auto&& self, [[maybe_unused]] size_t index) noexcept {
        assert(index == 0 && "[Error]: Index overflow");
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, 1, 1>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        if constexpr (v.isLValueMatrix())
            return LMatrixBlock<M1, 1, 1>(std::forward<M1>(v), self.fromRow, self.fromCol);
        else
            return RMatrixBlock<M1, 1, Dynamic>(std::forward<M1>(v), self.fromRow, self.fromCol, 1);
    }

    template<Matrix M>
    class LMatrixBlock<M, Dynamic, Dynamic> : public LValueMatrix<LMatrixBlock<M, Dynamic, Dynamic>> {
        using This = LMatrixBlock<M, Dynamic, Dynamic>;
        using Base = LValueMatrix<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }

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
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
    };

    template<Matrix M>
    LMatrixBlock<M, Dynamic, Dynamic>::LMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        Base::checkBlock(mat, fromRow, rowCount, fromCol, colCount);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, 1, Dynamic>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::topRows(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::bottomRows(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, 0, self.getRow(), fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::leftCols(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, self.getRow(), 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::rightCols(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, 0, self.getRow(), from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        Base::checkBlock(self, 0, toRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::topLeftCorner(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        Base::checkBlock(self, 0, toRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.rowCount);
        assert(col < self.colCount);
        return self.mat.data_ptr(row + self.fromRow, col + self.fromCol);
    }

    template<Matrix M>
    auto LMatrixBlock<M, Dynamic, Dynamic>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        if constexpr (v.isLValueMatrix())
            return LMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol, self.colCount);
        else
            return RMatrixBlock<M1, Dynamic, Dynamic>(std::forward<M1>(v), self.fromRow, self.rowCount, self.fromCol, self.colCount);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<LMatrixBlock<M, Row, Col>> {
        static_assert(std::remove_cvref_t<M>::isLValueMatrix());
        using M1 = std::remove_cvref<M>::type;
    public:
        using ScalarType = M1::ScalarType;
        constexpr static int Major = M1::Major;
    };
}
