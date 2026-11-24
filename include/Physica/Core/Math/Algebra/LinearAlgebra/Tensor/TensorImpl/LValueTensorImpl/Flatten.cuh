/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.cuh"
#include "../LValueTensor.cuh"
#include "Flatten.h"

namespace Physica {
    template<Tensor T> requires(std::remove_cvref_t<T>::isLValueTensor() && !std::remove_cvref_t<T>::isCompact())
    class device_obj<Flatten<T>> : public device_obj<LValueVector<Flatten<T>>> {
        using host_obj = Flatten<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
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
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return tensor.getDerived().getSize(); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;
    };

    template<Tensor T> requires(std::remove_cvref_t<T>::isLValueTensor() && !std::remove_cvref_t<T>::isCompact())
    __host__ __device__ void device_obj<Flatten<T>>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Tensor T> requires(std::remove_cvref_t<T>::isLValueTensor() && !std::remove_cvref_t<T>::isCompact())
    __host__ __device__ auto device_obj<Flatten<T>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.tensor.getDerived().data_ptr(self.tensor.getDerived().toIndexND(index));
    }
}
