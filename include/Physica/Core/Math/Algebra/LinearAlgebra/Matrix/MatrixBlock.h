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
namespace Physica::Core {
    template<class MatrixType>
    class MatrixBlock {
        MatrixType& mat;
        size_t fromRow;
        size_t toRow;
        size_t fromCol;
        size_t toCol;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t toRow_, size_t fromCol_, size_t toCol_);
        MatrixBlock(const MatrixBlock&) = delete;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        MatrixBlock& operator=(const MatrixBlock&) = delete;
        MatrixBlock& operator=(MatrixBlock&&) noexcept = delete;
        T& operator()(size_t row, size_t col) { assert((row + fromRow) < toRow); return mat(row + fromRow, col + fromCol); }
        const T& operator()(size_t row, size_t col) const { assert((col + fromCol) < toCol); return mat(row + fromRow, col + fromCol); }
    };

    template<class MatrixType>
    MatrixBlock<class MatrixType>::MatrixBlock(MatrixType& vec_, size_t fromRow_, size_t toRow_, size_t fromCol_, size_t toCol_)
            : mat(mat_)
            , fromRow(fromRow_)
            , toRow(toRow_)
            , fromCol(fromCol_)
            , toCol(toCol_) {
        assert(fromRow < toRow);
        assert(fromCol < toCol);
        assert(toRow <= mat.getRow());
        assert(toCol <= mat.getColumn());
    }
}