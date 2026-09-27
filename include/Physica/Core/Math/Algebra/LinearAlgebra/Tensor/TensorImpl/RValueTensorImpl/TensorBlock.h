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

namespace Physica {
    /**
     * Reference a hyper-rectangular part of the given tensor
     */
    template<class X>
    class TensorBlock : public RValueTensor<TensorBlock<X>> {
        using This = TensorBlock<X>;
        using Base = RValueTensor<This>;
    public:
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
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
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::resize;
        void resize(IndexType size);
        [[nodiscard]] T calc(const IndexType& index) const;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return shape; }
    };

    template<class X>
    TensorBlock<X>::TensorBlock(X&& tensor_, IndexType from_, IndexType count)
            : tensor(std::forward<X>(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor.dim(i));
            assert(from[i] + shape[i] <= tensor.dim(i));
        }
    }

    template<class X>
    void TensorBlock<X>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<class X>
    auto TensorBlock<X>::calc(const IndexType& index) const -> T {
        IndexType global{};
        for (auto&& [g, f, i] : zip(global, from, index))
            g = f + i;
        return tensor.calc(global);
    }
}

namespace Physica {
    template<class X>
    class Traits<TensorBlock<X>> {
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
        constexpr static int NDim = std::remove_cvref_t<X>::NDim;
    };
}
