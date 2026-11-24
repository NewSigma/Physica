/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.cuh"
#include "../RValueTensor.cuh"
#include "TensorFiber.h"

namespace Physica {
    template<class X, int Dim>
    class device_obj<TensorFiber<X, Dim>> : public device_obj<RValueVector<TensorFiber<X, Dim>>> {
        using host_obj = TensorFiber<X, Dim>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj_t<X>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<X>>> tensor;
        IndexType index;
    public:
        __host__ __device__ device_obj(Ref tensor, IndexVar auto... indices);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return tensor.getDerived().dim(Dim); }
    };

    template<class X, int Dim>
    __host__ __device__ device_obj<TensorFiber<X, Dim>>::device_obj(Ref tensor_, IndexVar auto... indices) : tensor(asStruct(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(static_cast<size_t>(indices) < tensor_.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<class X, int Dim>
    __device__ auto device_obj<TensorFiber<X, Dim>>::calc(size_t i, [[maybe_unused]] instanceof_x<ThreadBlock> auto block) const -> T {
        assert(i < getLength());
        auto idx = index;
        idx[Dim] = i;
        return tensor.getDerived().calc(idx);
    }
}

namespace Physica {
    template<class X, int Dim>
    class Traits<device_obj<TensorFiber<X, Dim>>> : public Traits<TensorFiber<X, Dim>> {};
}
