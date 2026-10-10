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

#include <functional>
#include <numeric>
#include "Physica/Macro.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    template<size_t Length = Dynamic>
    class Index : public Array<size_t, Length> {
        using This = Index<Length>;
        using Base = Array<size_t, Length>;
    public:
        using Base::Base;
        Index(const Base& base) noexcept : Base(base) {}
        Index(const This&) = default;
        Index(This&&) noexcept = default;
        ~Index() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Static members */
        [[nodiscard]] __host__ __device__ static size_t toIndex1D(const This& __restrict shape, const This& __restrict indices) noexcept;
        [[nodiscard]] __host__ __device__ static This toIndexND(const This& shape, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ static size_t toSize(const This& shape) noexcept;
    };

    template<size_t Length>
    __host__ __device__ size_t Index<Length>::toIndex1D(const This& __restrict shape, const This& __restrict indices) noexcept {
        size_t index = 0;
        size_t stride = 1;
        for (int i = static_cast<int>(shape.getLength()) - 1; i >= 0; --i) {
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            index += indices[i] * stride;
            stride *= shape[i];
        }
        return index;
    }

    template<size_t Length>
    __host__ __device__ auto Index<Length>::toIndexND(const This& shape, size_t index) noexcept -> This {
        static_assert(Length != Dynamic, "[Error]: Dynamic index can not decode a 1D index on device");
        const int Dim = shape.getLength();
        This indices(Dim);
        size_t stride = 1;
        for (int i = Dim - 1; i >= 0; --i) {
            indices[i] = (index / stride) % shape[i];
            stride *= shape[i];
        }
        assert(index < stride && "[Error]: Index out of range");
        return indices;
    }

    template<size_t Length>
    __host__ __device__ size_t Index<Length>::toSize(const This& shape) noexcept {
        return std::reduce(shape.begin(), shape.end(), size_t{1}, std::multiplies<>{});
    }

    template<size_t Length>
    void forND(const Index<Length>& shape, std::invocable<Index<Length>> auto func) {
        Index<Length> index(shape.getLength(), 0);
        [&](this auto&& self, size_t dim) -> void {
            if (dim == shape.getLength())
                func(index);
            else {
                for (size_t i = 0; i < shape[dim]; ++i) {
                    index[dim] = i;
                    self(dim + 1);
                }
            }
        }(0);
    }

    using Index1D = Index<1>;
    using Index2D = Index<2>;
    using Index3D = Index<3>;
    using Index4D = Index<4>;
    using Index5D = Index<5>;
    using IndexND = Index<>;
}
