/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../RValueTensor.cuh"
#include "TensorBlock.h"

namespace Physica {
    template<class X>
    class device_obj<TensorBlock<X>> : public device_obj<RValueTensor<TensorBlock<X>>> {
        using host_obj = TensorBlock<X>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
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
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(IndexType size);
        [[nodiscard]] __device__ T calc(const IndexType& index) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t dim(int index) const noexcept { return shape[index]; }
        [[nodiscard]] __host__ __device__ constexpr IndexType getShape() const noexcept { return shape; }
    };

    template<class X>
    __host__ __device__ device_obj<TensorBlock<X>>::device_obj(Ref tensor_, IndexType from_, IndexType count)
            : tensor(asStruct(tensor_))
            , from(std::move(from_))
            , shape(std::move(count)) {
        for (int i = 0; i < Base::NDim; ++i) {
            assert(from[i] < tensor_.dim(i));
            assert(from[i] + shape[i] <= tensor_.dim(i));
        }
    }

    template<class X>
    __host__ __device__ void device_obj<TensorBlock<X>>::resize([[maybe_unused]] IndexType size) {
        assert(size == shape && "[Error]: Resize part of a grid is not allowed");
    }

    template<class X>
    __device__ auto device_obj<TensorBlock<X>>::calc(const IndexType& index) const -> T {
        IndexType global{};
        for (int i = 0; i < Base::NDim; ++i)
            global[i] = from[i] + index[i];
        return tensor.getDerived().calc(global);
    }
}

namespace Physica {
    template<class X>
    class Traits<device_obj<TensorBlock<X>>> : public Traits<TensorBlock<X>> {};
}
