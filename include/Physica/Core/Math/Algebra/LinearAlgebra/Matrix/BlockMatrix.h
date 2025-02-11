/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Utils/Unreachable.h"
#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Matrix T>
    class BlockMatrix : public RValueMatrix<BlockMatrix<T>> {
        using This = BlockMatrix<T>;
        using Base = RValueMatrix<This>;
        using BlockArray = Array<T>;
    public:
        using typename Base::ScalarType;
    private:
        BlockArray blocks;
        Array<size_t> indexEnds;
    public:
        BlockMatrix(BlockArray blocks_);
        BlockMatrix(const This&) = default;
        BlockMatrix(This&&) noexcept = default;
        ~BlockMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const BlockArray& getBlocks() const noexcept { return blocks; }
        [[nodiscard]] size_t getNumBlock() const noexcept { return blocks.getLength(); }
        [[nodiscard]] const Array<size_t>& getIndexEnds() const noexcept { return indexEnds; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return indexEnds[getNumBlock() - 1]; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getRow(); }
    private:
        void updateEnds();
        size_t findBlock(size_t globalIndex) const;
    };

    template<Matrix T>
    BlockMatrix<T>::BlockMatrix(BlockArray blocks_) : blocks(std::move(blocks_)) {
        assert(getNumBlock() != 0 && "[Error]: Empty blocks does nothing");
        updateEnds();
    }

    template<Matrix T>
    BlockMatrix<T>::ScalarType BlockMatrix<T>::calc(size_t row, size_t col) const {
        const size_t indexR = findBlock(row);
        const size_t indexC = findBlock(col);
        if (indexR != indexC)
            return ScalarType(0);
        const size_t shift = indexR == 0 ? static_cast<size_t>(0) : indexEnds[indexR - 1];
        return blocks[indexR].calc(row - shift, col - shift);
    }

    template<Matrix T>
    void BlockMatrix<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        blocks.swap(obj.blocks);
        indexEnds.swap(obj.indexEnds);
    }

    template<Matrix T>
    void BlockMatrix<T>::updateEnds() {
        const size_t numBlock = getNumBlock();
        indexEnds.resize(numBlock);
        size_t end = 0;
        for (size_t i = 0; i < numBlock; ++i) {
            end += blocks[i].getRow();
            indexEnds[i] = end;
        }
    }

    template<Matrix T>
    size_t BlockMatrix<T>::findBlock(size_t globalIndex) const {
        assert(globalIndex < getRow() && "[Error]: Index out of range");
        for (size_t i = 0; i < getNumBlock(); ++i) {
            const bool inTheBlock = globalIndex < indexEnds[i];
            if (inTheBlock)
                return i;
        }
        unreachable();
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<BlockMatrix<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#include "BlockMatrixImpl/MatrixVectorProduct.h"
