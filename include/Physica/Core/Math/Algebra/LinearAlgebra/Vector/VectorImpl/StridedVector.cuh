/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "LValueVector.cuh"
#include "StridedVector.h"

namespace Physica {
    template<class Derived>
    class device_obj<StridedVector<Derived>> : public device_obj<LValueVector<Derived>> {
        using host_obj = StridedVector<Derived>;
        using Base = device_obj<LValueVector<Derived>>;
        using This = device_obj<host_obj>;
    protected:
        using typename Base::T;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        device_obj<Derived>& operator=(Scalar auto x);
        using Base::operator=;
        /* Operations */
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(this auto&&, size_t to = Length) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(this auto&&, size_t from, size_t to) noexcept;

        void zeros();
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedVector<Derived>>::getStride() const noexcept {
        if constexpr (Derived::getStrideAtCompile() != Dynamic)
            return Derived::getStrideAtCompile();
        else
            return Base::getDerived().getStride();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.data_handle() + index * self.getStride();
    }
}

#include "StridedVectorImpl/StridedVectorImpl.cuh"
#include "StridedVectorImpl/StridedVectorBlock.cuh"
#include "StridedVectorImpl/VectorConvert/RealVector.cuh"
#include "StridedVectorImpl/VectorConvert/ImagVector.cuh"
