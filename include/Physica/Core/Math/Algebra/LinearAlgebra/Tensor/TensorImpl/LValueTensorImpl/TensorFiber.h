/*
 * Copyright 2025 Weibo He.
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
    template<Tensor X>
    class TensorFiber : public LValueVector<TensorFiber<X>> {
        using This = TensorFiber<X>;
        using Base = LValueVector<TensorFiber<X>>;
        using IndexType = X::IndexType;
        constexpr static int NDim = X::NDim;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        X& tensor;
        IndexType index;
        int dim;
    public:
        TensorFiber(X& tensor, int dim, IndexType index);
        TensorFiber(const This&) = default;
        TensorFiber(This&&) noexcept = default;
        ~TensorFiber() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t length);
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return tensor.dim(dim); }
        [[nodiscard]] PtrTy data_ptr(size_t i) noexcept;
        using Base::data_ptr;
    };

    template<Tensor X>
    TensorFiber<X>::TensorFiber(X& tensor, int dim, IndexType index)
            : tensor(tensor), index(index), dim(dim) {
        assert(dim < NDim);
        for (int i = 0; i < NDim; ++i) {
            if (i == dim)
                continue;
            assert(index[i] < tensor.dim(i));
        }
    }

    template<Tensor X>
    void TensorFiber<X>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Tensor X>
    auto TensorFiber<X>::data_ptr(size_t i) noexcept -> PtrTy {
        auto idx = index;
        idx[dim] = i;
        return tensor.data_ptr(idx);
    }
}

namespace Physica {
    template<Tensor X>
    class Traits<TensorFiber<X>> {
    public:
        using ScalarType = X::ScalarType;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
