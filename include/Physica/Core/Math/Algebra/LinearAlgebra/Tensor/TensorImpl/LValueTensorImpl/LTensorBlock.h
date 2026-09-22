/*
 * Copyright 2023-2026 Weibo He.
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
    template<class Derived> class LValueTensor;
    template<class TensorType> class LTensorBlock;

    template<Tensor X>
    class LTensorBlock<X> : public LValueTensor<LTensorBlock<X>> {
        using This = LTensorBlock<X>;
        using Base = LValueTensor<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::IndexType;
        using Base::NDim;
    private:
        decay_rvalue_t<X> grid;
        IndexType from;
        IndexType shape;
    public:
        LTensorBlock(X&& grid_, IndexType from, IndexType count);
        LTensorBlock(const LTensorBlock&) = delete;
        LTensorBlock(LTensorBlock&&) noexcept = delete;
        ~LTensorBlock() = default;
        /* Operators */
        using Base::operator=;
        LTensorBlock& operator=(const LTensorBlock& b) { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        LTensorBlock& operator=(LTensorBlock&& b) noexcept { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        /* Operations */
        using Base::resize;
        void resize(IndexType size);

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return shape; }
        [[nodiscard]] auto data_ptr(this auto&&, const IndexType& index) noexcept;
    };

    template<Tensor X>
    LTensorBlock<X>::LTensorBlock(X&& grid_, IndexType from, IndexType count)
            : grid(std::forward<X>(grid_))
            , from(std::move(from))
            , shape(std::move(count)) {
        for (int i = 0; i < NDim; ++i) {
            assert(from[i] < grid.dim(i));
            assert(from[i] + count[i] <= grid.dim(i));
        }
    }

    template<Tensor X>
    void LTensorBlock<X>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<Tensor X>
    auto LTensorBlock<X>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        IndexType global{};
        for (auto&& [g, f, i] : zip(global, self.from, index))
            g = f + i;
        return self.grid.data_ptr(global);
    }

    template<Tensor X>
    auto LTensorBlock<X>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), X>(self.grid).values();
        using X1 = decltype(v);
        return LTensorBlock<X1>(std::forward<X1>(v), self.from, self.shape);
    }
}

namespace Physica {
    template<Tensor X>
    class Traits<LTensorBlock<X>> {
        static_assert(std::remove_cvref_t<X>::isLValueTensor());
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
        constexpr static int NDim = std::remove_cvref_t<X>::NDim;
    };
}
