/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "CompactVector.h"
#include "StridedVector.cuh"
#include "CompactVectorImpl/CompactVectorBlock.cuh"

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
        device_obj<Derived>& operator=(Scalar auto x);
        using Base::operator=;
        using Base::operator+=;
        /* Operations */
        template<int Size>
        [[nodiscard]] __device__ auto packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ auto packet(size_t index, size_t count) const noexcept;
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

        void zeros();
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();

        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data() noexcept;
        [[nodiscard]] __host__ __device__ auto data() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastPacket() noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "CompactVectorImpl/CompactVectorImpl.cuh"
#include "CompactVectorImpl/CompactVectorBlock.cuh"
#include "CompactVectorImpl/VectorConvert/RealVector.cuh"
#include "CompactVectorImpl/VectorConvert/ImagVector.cuh"
