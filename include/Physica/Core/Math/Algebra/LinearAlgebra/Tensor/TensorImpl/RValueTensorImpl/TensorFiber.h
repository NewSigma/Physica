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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
    template<class X, int Dim>
    class TensorFiber : public RValueVector<TensorFiber<X, Dim>> {
        using This = TensorFiber<X, Dim>;
        using Base = RValueVector<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
        IndexType index;
    public:
        TensorFiber(X&& tensor, IndexVar auto... indices);
        TensorFiber(const This&) = default;
        TensorFiber(This&&) noexcept = default;
        ~TensorFiber() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t i) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return tensor.dim(Dim); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<class X, int Dim>
    TensorFiber<X, Dim>::TensorFiber(X&& tensor_, IndexVar auto... indices) : tensor(std::forward<X>(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<class X, int Dim>
    auto TensorFiber<X, Dim>::calc(size_t i) const -> T {
        assert(i < getLength());
        auto idx = index;
        idx[Dim] = i;
        return tensor.calc(idx);
    }
}

namespace Physica {
    template<class X, int Dim>
    class Traits<TensorFiber<X, Dim>> {
    public:
        using ScalarType = std::remove_cvref_t<X>::ScalarType;
    };
}
