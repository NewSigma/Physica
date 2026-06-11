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
    private:
        decay_rvalue_t<X> grid;
        Index3D from;
        Index3D count;
    public:
        LTensorBlock(X&& grid_, Index3D from, Index3D count);
        LTensorBlock(const LTensorBlock&) = delete;
        LTensorBlock(LTensorBlock&&) noexcept = delete;
        ~LTensorBlock() = default;
        /* Operators */
        using Base::operator=;
        LTensorBlock& operator=(const LTensorBlock& b) { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        LTensorBlock& operator=(LTensorBlock&& b) noexcept { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        /* Operations */
        void resize([[maybe_unused]] Index3D size) { assert(size == count && "[Error]: Resize part of a grid is not allowed"); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return count[0]; }
        [[nodiscard]] size_t getDimY() const noexcept { return count[1]; }
        [[nodiscard]] size_t getDimZ() const noexcept { return count[2]; }
        [[nodiscard]] auto data_ptr(this auto&&, Index3D index) noexcept;
    };

    template<Tensor X>
    LTensorBlock<X>::LTensorBlock(X&& grid_, Index3D from, Index3D count)
            : grid(std::forward<X>(grid_))
            , from(from)
            , count(count) {
        for (int i = 0; i < 3; ++i) {
            assert(from[i] < grid.getDim()[i]);
            assert(from[i] + count[i] <= grid.getDim()[i]);
        }
    }

    template<Tensor X>
    auto LTensorBlock<X>::data_ptr(this auto&& self, Index3D index) noexcept {
        return self.grid.data_ptr({self.from[0] + index[0], self.from[1] + index[1], self.from[2] + index[2]});
    }

    template<Tensor X>
    auto LTensorBlock<X>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), X>(self.grid).values();
        using X1 = decltype(v);
        return LTensorBlock<X1>(std::forward<X1>(v), self.from, self.count);
    }
}

namespace Physica {
    template<Tensor X>
    class Traits<LTensorBlock<X>> {
        static_assert(std::remove_cvref_t<X>::isLValueTensor());
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
    };
}
