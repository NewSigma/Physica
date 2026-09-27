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

#include "../LValueTensor.h"

namespace Physica {
    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    class TensorBlock<X> : public LValueTensor<TensorBlock<X>> {
        using This = TensorBlock<X>;
        using Base = LValueTensor<This>;
    public:
        using typename Base::ScalarType;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    private:
        decay_rvalue_t<X> tensor;
        IndexType from;
        IndexType shape;
    public:
        TensorBlock(X&& tensor, IndexType from, IndexType count);
        TensorBlock(const This&) = delete;
        TensorBlock(This&&) noexcept = delete;
        ~TensorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& b);
        This& operator=(This&& b) noexcept;
        /* Operations */
        using Base::resize;
        void resize(IndexType size);
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return shape; }
        [[nodiscard]] auto data_ptr(this auto&&, const IndexType& index) noexcept;
    };

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    TensorBlock<X>::TensorBlock(X&& tensor_, IndexType from_, IndexType count)
            : tensor(std::forward<X>(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor.dim(i));
            assert(from[i] + shape[i] <= tensor.dim(i));
        }
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::operator=(const This& b) -> This& {
        Base::operator=(static_cast<const Base::Base&>(b));
        return *this;
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::operator=(This&& b) noexcept -> This& {
        Base::operator=(static_cast<const Base::Base&>(b));
        return *this;
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    void TensorBlock<X>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), X>(self.tensor).values();
        using X1 = decltype(v);
        return TensorBlock<X1>(std::forward<X1>(v), self.from, self.shape);
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    auto TensorBlock<X>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        IndexType global{};
        for (auto&& [g, f, i] : zip(global, self.from, index))
            g = f + i;
        return self.tensor.data_ptr(global);
    }
}
