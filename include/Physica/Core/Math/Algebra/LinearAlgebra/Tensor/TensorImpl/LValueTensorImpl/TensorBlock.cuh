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

#include "../LValueTensor.cuh"
#include "TensorBlock.h"

namespace Physica {
    template<class X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    class device_obj<TensorBlock<X>> : public device_obj<LValueTensor<TensorBlock<X>>> {
        using host_obj = TensorBlock<X>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueTensor<host_obj>>;
        using Ref = add_device_obj_t<X>;
    public:
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<X>>> tensor;
        IndexType from;
        IndexType shape;
    public:
        __host__ __device__ device_obj(Ref tensor, IndexType from, IndexType count);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(IndexType size);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return shape[index]; }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return shape; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
    };

    template<class X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    __host__ __device__ device_obj<TensorBlock<X>>::device_obj(Ref tensor_, IndexType from_, IndexType count)
            : tensor(asStruct(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor_.dim(i));
            assert(from[i] + shape[i] <= tensor_.dim(i));
        }
    }

    template<class X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    __host__ __device__ void device_obj<TensorBlock<X>>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<class X> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isCompact())
    __host__ __device__ auto device_obj<TensorBlock<X>>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        IndexType global{};
        for (int i = 0; i < Base::NDim; ++i)
            global[i] = self.from[i] + index[i];
        return self.tensor.getDerived().data_ptr(global);
    }
}
