/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "CompactVector.h"
#include "CompactVectorImpl/CompactVectorBlock.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "StridedVector.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<CompactVector<Derived>> : public device_obj<StridedVector<Derived>> {
        using host_obj = CompactVector<Derived>;
        using Base = device_obj<StridedVector<Derived>>;
        using This = device_obj<host_obj>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        using DataSetType = host_obj::DataSetType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        using Base::operator+=;
        /* Operations */
        template<int Size>
        [[nodiscard]] __device__ device_obj<SIMD<T, Size>> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ device_obj<SIMD<T, Size>> packet(size_t index, size_t count) const noexcept;
        __device__ void writePacket(Packet auto packet, size_t index) noexcept;
        __device__ void writePacket(Packet auto packet, size_t index, size_t count) noexcept;
        void reverse(const auto& grad) const noexcept;

        template<Vector V> void toHost(CompactVector<V>& obj) const;
        template<Vector V> void toHostAsync(CompactVector<V>& obj) const;

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(this auto&&, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(this auto&&, size_t from, size_t to) noexcept;

        void read(const auto& obj);
        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "CompactVectorImpl/CompactVectorImpl.cuh"
#include "CompactVectorImpl/CompactVectorBlock.cuh"
