/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <functional>
#include <numeric>
#include "Physica/Macro.h"
#include "Physica/Core/Utils/Builtin.h"
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
        /* Operations */
        using Base::zeros;
        /* Static members */
        [[nodiscard]] __host__ __device__ static size_t toIndex1D(const This& __restrict shape, const This& __restrict indices) noexcept;
        [[nodiscard]] __host__ __device__ static This toIndexND(const This& shape, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ static size_t toSize(const This& shape) noexcept;
        [[nodiscard]] __host__ __device__ constexpr static This zeros(size_t length) noexcept;
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
        return std::reduce(shape.begin(), shape.end(), size_t{1}, Physica::ckd_mul);
    }

    template<size_t Length>
    __host__ __device__ constexpr auto Index<Length>::zeros(size_t length) noexcept -> This {
        if constexpr (Length == Dynamic)
            return This(length, size_t{0});
        else {
            assert(length == Length && "[Error]: Length mismatch");
            return This{};
        }
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
