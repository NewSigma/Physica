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
    template<Tensor X, int DimR, int DimC>
    class TensorSlice : public LValueMatrix<TensorSlice<X, DimR, DimC>> {
        using This = TensorSlice<X, DimR, DimC>;
        using Base = LValueMatrix<TensorSlice<X, DimR, DimC>>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(DimR != DimC, "[Error]: DimR and DimC must be different");
        static_assert(DimR < std::remove_cvref_t<X>::ndim());
        static_assert(DimC < std::remove_cvref_t<X>::ndim());
    private:
        decay_rvalue_t<X> tensor;
        IndexType index;
    public:
        TensorSlice(X&& tensor, IndexVar auto... indices);
        TensorSlice(const This&) = default;
        TensorSlice(This&&) noexcept = default;
        ~TensorSlice() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tensor.dim(DimR); }
        [[nodiscard]] size_t getCol() const noexcept { return tensor.dim(DimC); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
    };

    template<Tensor X, int DimR, int DimC>
    TensorSlice<X, DimR, DimC>::TensorSlice(X&& tensor, IndexVar auto... indices)
            : tensor(std::forward<X>(tensor)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<Tensor X, int DimR, int DimC>
    void TensorSlice<X, DimR, DimC>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow() && col == getCol());
    }

    template<Tensor X, int DimR, int DimC>
    auto TensorSlice<X, DimR, DimC>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        auto idx = self.index;
        idx[DimR] = row;
        idx[DimC] = col;
        return self.tensor.data_ptr(idx);
    }
}

namespace Physica {
    template<Tensor X, int DimR, int DimC>
    class Traits<TensorSlice<X, DimR, DimC>> {
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
    };
}
