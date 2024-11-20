/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    template<class Derived>
    class device_obj<ContinuousVector<Derived>> : public device_obj<LValueVector<Derived>> {
        using Base = device_obj<LValueVector<Derived>>;
        using This = device_obj<ContinuousVector<Derived>>;
        template<size_t Length> using BlockType = device_obj<ContinuousVectorBlock<Derived, Length>>;
    public:
        using host_obj = ContinuousVector<Derived>;
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        ~device_obj() = default;
        /* Operators */
        inline This& operator=(const This& obj);
        inline This& operator=(This&& obj);
        using Base::operator=;
        /* Operations */
        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ inline AnyPacket packet(size_t index) const;
        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ inline AnyPacket packetPartial(size_t index, size_t count) const;
        template<class AnyPacket> __device__ inline void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> __device__ inline void writePacketPartial(size_t index, size_t count, const AnyPacket packet);

        void resize(size_t length) { Base::getDerived().resize(length); }
        template<Vector V> void toHost(ContinuousVector<V>& obj) const;
        template<Vector V> void toHostAsync(ContinuousVector<V>& obj) const;

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> segment(size_t from, size_t to) const;
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
