/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Utils/Builtin.h"
#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class BlockMatrix : public RValueMatrix<BlockMatrix<M>> {
        using This = BlockMatrix<M>;
        using Base = RValueMatrix<This>;
        using BlockArray = Array<M>;
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
        [[nodiscard]] const auto& getBlocks() const noexcept { return blocks; }
        [[nodiscard]] size_t getNumBlock() const noexcept { return blocks.getLength(); }
        [[nodiscard]] const auto& getIndexEnds() const noexcept { return indexEnds; }
        [[nodiscard]] size_t getOrder() const noexcept { return indexEnds[getNumBlock() - 1]; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return true; }
    private:
        void updateEnds();
        [[nodiscard]] size_t findBlock(size_t globalIndex) const;
    };

    template<Matrix M>
    BlockMatrix<M>::BlockMatrix(BlockArray blocks_) : blocks(std::move(blocks_)) {
        assert(getNumBlock() != 0 && "[Error]: Empty blocks does nothing");
        updateEnds();
    }

    template<Matrix M>
    BlockMatrix<M>::ScalarType BlockMatrix<M>::calc(size_t row, size_t col) const {
        const size_t indexR = findBlock(row);
        const size_t indexC = findBlock(col);
        if (indexR != indexC)
            return ScalarType(0);
        const size_t shift = indexR == 0 ? static_cast<size_t>(0) : indexEnds[indexR - 1];
        return blocks[indexR].calc(row - shift, col - shift);
    }

    template<Matrix M>
    void BlockMatrix<M>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        blocks.swap(obj.blocks);
        indexEnds.swap(obj.indexEnds);
    }

    template<Matrix M>
    void BlockMatrix<M>::updateEnds() {
        const size_t numBlock = getNumBlock();
        indexEnds.resize(numBlock);
        size_t end = 0;
        for (size_t i = 0; i < numBlock; ++i) {
            end += blocks[i].getRow();
            indexEnds[i] = end;
        }
    }

    template<Matrix M>
    size_t BlockMatrix<M>::findBlock(size_t globalIndex) const {
        assert(globalIndex < getOrder() && "[Error]: Index out of range");
        for (size_t i = 0; i < getNumBlock(); ++i) {
            const bool inTheBlock = globalIndex < indexEnds[i];
            if (inTheBlock)
                return i;
        }
        unreachable();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<BlockMatrix<M>> {
    public:
        using ScalarType = M::ScalarType;
        constexpr static int Prec = M::Prec;
    };
}

#include "BlockMatrixImpl/GEMV.h"
