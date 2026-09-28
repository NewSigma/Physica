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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.cuh"
#include "../LValueTensor.cuh"
#include "TensorFiber.h"

namespace Physica {
    template<class X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    class device_obj<TensorFiber<X, Dim>> : public device_obj<LValueVector<TensorFiber<X, Dim>>> {
        using host_obj = TensorFiber<X, Dim>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using Ref = add_device_obj_t<X>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(Dim < std::remove_cvref_t<X>::ndim());
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
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return tensor.getDerived().dim(Dim); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t i) noexcept;
    };

    template<class X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    __host__ __device__ device_obj<TensorFiber<X, Dim>>::device_obj(Ref tensor_, IndexVar auto... indices) : tensor(asStruct(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor_.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<class X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    __host__ __device__ void device_obj<TensorFiber<X, Dim>>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<class X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    __host__ __device__ auto device_obj<TensorFiber<X, Dim>>::data_ptr(this auto&& self, size_t i) noexcept {
        assert(i < self.getLength());
        auto idx = self.index;
        idx[Dim] = i;
        return self.tensor.getDerived().data_ptr(idx);
    }
}
