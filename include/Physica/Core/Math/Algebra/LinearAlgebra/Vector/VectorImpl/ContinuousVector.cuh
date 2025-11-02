/*
 * Copyright 2023-2025 Weibo He.
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
        template<size_t Length>
        using BlockType = device_obj<ContinuousVectorBlock<Derived, Length>>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
        using Base::MaxThreadPerBlock;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using DataSetType = host_obj::DataSetType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) = delete;
        using Base::operator=;
        using Base::operator+=;
        /* Operations */
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const;
        __device__ void writePacket(size_t index, const Packet auto packet);
        __device__ void writePacketPartial(size_t index, size_t count, const Packet auto packet);
        void reverse(const auto& grad) const noexcept requires(isReverseDiff);

        template<Vector V> void toHost(ContinuousVector<V>& obj) const;
        template<Vector V> void toHostAsync(ContinuousVector<V>& obj) const;

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto segment(size_t from, size_t to) const noexcept;

        void zeros();

        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data() { return Base::data_ptr(0); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data() const { return Base::data_ptr(0); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.cuh"
#include "ContinuousVectorImpl/ContinuousVectorBlock.cuh"
