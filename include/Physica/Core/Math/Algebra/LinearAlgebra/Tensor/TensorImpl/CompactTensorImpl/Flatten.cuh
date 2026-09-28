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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.cuh"
#include "../CompactTensor.cuh"
#include "Flatten.h"

namespace Physica {
    template<Tensor T> requires(std::remove_cvref_t<T>::isCompact())
    class device_obj<Flatten<T>> : public device_obj<CompactVector<Flatten<T>>> {
        using host_obj = Flatten<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;
        using Ref = add_device_obj_t<T>;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<T>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
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
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return tensor.getDerived().getSize(); }
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
    };

    template<Tensor T> requires(std::remove_cvref_t<T>::isCompact())
    __host__ __device__ void device_obj<Flatten<T>>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Tensor T> requires(std::remove_cvref_t<T>::isCompact())
    __host__ __device__ auto device_obj<Flatten<T>>::data_handle(this auto&& self) noexcept {
        return self.tensor.getDerived().data_handle();
    }
}
