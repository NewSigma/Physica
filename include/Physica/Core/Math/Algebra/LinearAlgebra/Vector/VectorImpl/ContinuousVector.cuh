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
#include "ContinuousVector.h"
#include "LValueVector.cuh"
#include "ContinuousVectorImpl/ContinuousVectorBlock.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<ContinuousVector<Derived>> : public device_obj<LValueVector<Derived>> {
        using host_obj = ContinuousVector<Derived>;
        using Base = device_obj<LValueVector<Derived>>;
        using This = device_obj<host_obj>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
        using Base::MaxThreadsPerBlock;
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
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const;
        __device__ void writePacket(size_t index, Packet auto packet);
        __device__ void writePacketPartial(size_t index, size_t count, Packet auto packet);
        void reverse(const auto& grad) const noexcept;

        template<Vector V> void toHost(ContinuousVector<V>& obj) const;
        template<Vector V> void toHostAsync(ContinuousVector<V>& obj) const;

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
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.cuh"
#include "ContinuousVectorImpl/ContinuousVectorBlock.cuh"
