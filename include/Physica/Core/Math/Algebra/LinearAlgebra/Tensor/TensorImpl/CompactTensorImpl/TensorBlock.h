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

#include "../CompactTensor.h"

namespace Physica {
    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    class TensorBlock<X> : public StridedTensor<TensorBlock<X>> {
        using This = TensorBlock<X>;
        using Base = StridedTensor<This>;
    public:
        using IndexType = std::remove_cvref_t<X>::IndexType;
    private:
        decay_rvalue_t<X> tensor;
        IndexType from;
        IndexType shape;
    public:
        TensorBlock(X&& tensor, IndexType from, IndexType count);
        TensorBlock(const This&) = default;
        TensorBlock(This&&) noexcept = default;
        ~TensorBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(IndexType size);
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return shape; }
        [[nodiscard]] auto getStrides() const noexcept;
        [[nodiscard]] auto data_handle(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
    };

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    TensorBlock<X>::TensorBlock(X&& tensor_, IndexType from_, IndexType count)
            : tensor(std::forward<X>(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor.dim(i));
            assert(from[i] + shape[i] <= tensor.dim(i));
        }
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    void TensorBlock<X>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), X>(self.tensor).values();
        using X1 = decltype(v);
        return TensorBlock<X1>(std::forward<X1>(v), self.from, self.shape);
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::getStrides() const noexcept {
        return tensor.getStrides();
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::data_handle(this auto&& self) noexcept {
        return self.tensor.data_ptr(self.from);
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval auto TensorBlock<X>::getStrideAtCompile() noexcept -> IndexType {
        return std::remove_cvref_t<X>::getStrideAtCompile();
    }
}
