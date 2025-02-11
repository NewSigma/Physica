/*
 * Copyright 2024 Weibo He.
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

#include "../ContinuousVector.cuh"

namespace Physica {
    template<class Derived>
    inline device_obj<ContinuousVector<Derived>>& device_obj<ContinuousVector<Derived>>::operator=(const This& obj) {
        Base::operator=(obj);
        return *this;
    }
    
    template<class Derived>
    inline device_obj<ContinuousVector<Derived>>& device_obj<ContinuousVector<Derived>>::operator=(This&& obj) {
        return *this = obj;
    }

    template<class Derived>
    template<Vector T>
    void device_obj<ContinuousVector<Derived>>::toHost(ContinuousVector<T>& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector T>
    void device_obj<ContinuousVector<Derived>>::toHostAsync(ContinuousVector<T>& obj) const {
        obj.resize(Base::getLength());
        check(cudaMemcpyAsync(obj.data(), data(), Base::getLength() * sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class Derived>
    template<class AnyPacket, Side Owner>
    __device__ inline AnyPacket device_obj<ContinuousVector<Derived>>::packet(size_t index) const {
        AnyPacket packet{};
        packet.load(Base::data_ptr(index));
        return packet;
    }

    template<class Derived>
    template<class AnyPacket, Side Owner>
    __device__ inline AnyPacket device_obj<ContinuousVector<Derived>>::packetPartial(size_t index, size_t count) const  {
        AnyPacket packet{};
        packet.load_partial(Base::data_ptr(index), count);
        return packet;
    }

    template<class Derived>
    template<class AnyPacket>
    __device__ inline void device_obj<ContinuousVector<Derived>>::writePacket(size_t index, const AnyPacket packet) {
        using LocalPacket = std::conditional<AnyPacket::size() == 1, ScalarType, SIMD<ScalarType, AnyPacket::size()>>::type;
        LocalPacket(packet).store(Base::data_ptr(index));
    }

    template<class Derived>
    template<class AnyPacket>
    __device__ inline void device_obj<ContinuousVector<Derived>>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        using LocalPacket = std::conditional<AnyPacket::size() == 1, ScalarType, SIMD<ScalarType, AnyPacket::size()>>::type;
        LocalPacket(packet).store_partial(Base::data_ptr(index), count);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline auto device_obj<ContinuousVector<Derived>>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<ContinuousVector<Derived>>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline auto device_obj<ContinuousVector<Derived>>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<ContinuousVector<Derived>>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline auto device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    template<Vector T>
    void ContinuousVector<Derived>::toDevice(device_obj<ContinuousVector<T>>& obj) const {
        toDeviceAsync(obj);
        if constexpr (!std::is_trivially_copy_constructible<T>::value)
            CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector T>
    void ContinuousVector<Derived>::toDeviceAsync(device_obj<ContinuousVector<T>>& obj) const {
        static_assert(std::is_same<ScalarType, typename T::ScalarType>::value,
                "[Error]: ScalarType inconsistent, additional buffer is necessary");
        const size_t length = Base::getLength();
        const size_t size = length * sizeof(ScalarType);
        obj.resize(length);
        if constexpr (std::is_trivially_copy_constructible<T>::value)
            memcpy(obj.data(), data(), size);
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyHostToDevice, CUDAContext::getInstance()));
    }
}
