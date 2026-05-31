/*
 * Copyright 2025-2026 Weibo He.
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

#include "../LValueTensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica {
    template<Tensor X>
    class TensorSlice : public LValueMatrix<TensorSlice<X>> {
        using This = TensorSlice<X>;
        using Base = LValueMatrix<TensorSlice<X>>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    private:
        LazyDestroy<X> tensor;
        IndexType index;
        int dimRow;
        int dimCol;
    public:
        TensorSlice(X&& tensor, int dimRow, int dimCol, IndexType index);
        TensorSlice(const This&) = default;
        TensorSlice(This&&) noexcept = default;
        ~TensorSlice() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tensor.dim(dimRow); }
        [[nodiscard]] size_t getCol() const noexcept { return tensor.dim(dimCol); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
    };

    template<Tensor X>
    TensorSlice<X>::TensorSlice(X&& tensor, int dimRow, int dimCol, IndexType index)
            : tensor(std::forward<X>(tensor)), index(index), dimRow(dimRow), dimCol(dimCol) {
        constexpr int NDim = tensor.ndim();
        assert(dimRow < NDim && dimCol < NDim);
        assert(dimRow != dimCol);
        for (int i = 0; i < NDim; ++i) {
            if (i == dimRow)
                continue;
            if (i == dimCol)
                continue;
            assert(index[i] < tensor.dim(i));
        }
    }

    template<Tensor X>
    void TensorSlice<X>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow() && col == getCol());
    }

    template<Tensor X>
    auto TensorSlice<X>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        auto idx = self.index;
        idx[self.dimRow] = row;
        idx[self.dimCol] = col;
        return self.tensor.data_ptr(idx);
    }
}

namespace Physica {
    template<Tensor X>
    class Traits<TensorSlice<X>> {
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
    };
}
