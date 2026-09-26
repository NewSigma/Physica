/*
 * Copyright 2026 Weibo He.
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

#include "../RValueTensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<class X, int DimR, int DimC>
    class TensorSlice : public RValueMatrix<TensorSlice<X, DimR, DimC>> {
        using This = TensorSlice<X, DimR, DimC>;
        using Base = RValueMatrix<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(DimR != DimC, "[Error]: DimR and DimC must be different");
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
        IndexType index;
    public:
        TensorSlice(X&& tensor, IndexVar auto... indices);
        TensorSlice(const This&) = default;
        TensorSlice(This&&) noexcept = default;
        ~TensorSlice() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tensor.dim(DimR); }
        [[nodiscard]] size_t getCol() const noexcept { return tensor.dim(DimC); }
        [[nodiscard]] size_t getOrder() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return DimR < DimC ? MatrixMajor::Row : MatrixMajor::Col; }
    };

    template<class X, int DimR, int DimC>
    TensorSlice<X, DimR, DimC>::TensorSlice(X&& tensor_, IndexVar auto... indices) : tensor(std::forward<X>(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<class X, int DimR, int DimC>
    auto TensorSlice<X, DimR, DimC>::calc(size_t row, size_t col) const -> T {
        assert(row < getRow());
        assert(col < getCol());
        auto idx = index;
        idx[DimR] = row;
        idx[DimC] = col;
        return tensor.calc(idx);
    }

    template<class X, int DimR, int DimC>
    size_t TensorSlice<X, DimR, DimC>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }
}

namespace Physica {
    template<class X, int DimR, int DimC>
    class Traits<TensorSlice<X, DimR, DimC>> {
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
    };
}
