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

#include "../CompactTensor.cuh"
#include "TensorBlock.h"

namespace Physica {
    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    class device_obj<TensorBlock<X>> : public device_obj<StridedTensor<TensorBlock<X>>> {
        using host_obj = TensorBlock<X>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedTensor<host_obj>>;
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
        [[nodiscard]] __device__ T calc(const IndexType& index) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return shape; }
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return shape[index]; }
        [[nodiscard]] __host__ __device__ auto getStrides() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
    };

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ device_obj<TensorBlock<X>>::device_obj(Ref tensor_, IndexType from_, IndexType count)
            : tensor(asStruct(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor_.dim(i));
            assert(from[i] + shape[i] <= tensor_.dim(i));
        }
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ void device_obj<TensorBlock<X>>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __device__ auto device_obj<TensorBlock<X>>::calc(const IndexType& index) const -> T {
        IndexType global{};
        for (int i = 0; i < Base::NDim; ++i)
            global[i] = from[i] + index[i];
        return tensor.getDerived().calc(global);
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ auto device_obj<TensorBlock<X>>::getStrides() const noexcept {
        return tensor.getDerived().getStrides();
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ auto device_obj<TensorBlock<X>>::data_handle(this auto&& self) noexcept {
        return self.tensor.getDerived().data_ptr(self.from);
    }

    template<Tensor X> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval auto device_obj<TensorBlock<X>>::getStrideAtCompile() noexcept -> IndexType {
        return std::remove_cvref_t<X>::getStrideAtCompile();
    }
}
