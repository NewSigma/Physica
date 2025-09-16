/*
 * Copyright 2024-2025 Weibo He.
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
    void device_obj<ContinuousVector<Derived>>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename T::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        if constexpr (Scalar<U>)
            Base::getConstCastDerived().grads() += grad;
        else {
            static_assert(CUDA<U>, "[Error]: Cannot pass host grad to device");
            if constexpr (Vector<U>) {
                assert(Base::getLength() == grad.getLength());
                Base::getConstCastDerived().grads() += grad;
            }
            else {
                static_assert(Matrix<U>, "[Error]: Unexpected type");
                assert(Base::getLength() == grad.getRow());
                reverse(grad.sum_cols());
            }
        }
    }

    template<class Derived>
    template<Vector V>
    void device_obj<ContinuousVector<Derived>>::toHost(ContinuousVector<V>& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector V>
    void device_obj<ContinuousVector<Derived>>::toHostAsync(ContinuousVector<V>& obj) const {
        static_assert(std::same_as<ScalarType, typename V::ScalarType>, "[Error]: Incompatible ScalarType");
        const size_t length = Base::getLength();
        const size_t size = length * sizeof(T);
        obj.resize(length);
        if constexpr (V::SizeAtCompile != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<T>) {
            Base::getDerived().values().toHostAsync(obj.getDerived().values());
            Base::getDerived().grads().toHostAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class Derived>
    template<Packet Pack>
    __device__ Pack device_obj<ContinuousVector<Derived>>::packet(size_t index) const {
        Pack packet{};
        if constexpr (isReverseDiff)
            packet.load(Base::data_ptr(index).value_ptr());
        else
            packet.load(Base::data_ptr(index));
        return packet;
    }

    template<class Derived>
    template<Packet Pack>
    __device__ Pack device_obj<ContinuousVector<Derived>>::packetPartial(size_t index, size_t count) const  {
        Pack packet{};
        if constexpr (isReverseDiff)
            packet.load_partial(Base::data_ptr(index).value_ptr(), count);
        else
            packet.load_partial(Base::data_ptr(index), count);
        return packet;
    }

    template<class Derived>
    template<Packet Pack>
    __device__ void device_obj<ContinuousVector<Derived>>::writePacket(size_t index, const Pack packet) {
        using T1 = std::conditional<isReverseDiff, Tv, T>::type;
        using LocalPacket = std::conditional<Pack::size() == 1, T1, SIMD<T1, Pack::size()>>::type;
        if constexpr (isReverseDiff)
            LocalPacket(packet).store(Base::data_ptr(index).value_ptr());
        else
            LocalPacket(packet).store(Base::data_ptr(index));
    }

    template<class Derived>
    template<Packet Pack>
    __device__ void device_obj<ContinuousVector<Derived>>::writePacketPartial(size_t index, size_t count, const Pack packet) {
        using T1 = std::conditional<isReverseDiff, Tv, T>::type;
        using LocalPacket = std::conditional<Pack::size() == 1, T1, SIMD<T1, Pack::size()>>::type;
        if constexpr (isReverseDiff)
            LocalPacket(packet).store_partial(Base::data_ptr(index).value_ptr(), count);
        else
            LocalPacket(packet).store_partial(Base::data_ptr(index), count);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<ContinuousVector<Derived>>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<ContinuousVector<Derived>>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<ContinuousVector<Derived>>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<ContinuousVector<Derived>>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    void device_obj<ContinuousVector<Derived>>::zeros() {
        if constexpr (Diffable<T>) {
            Base::getDerived().values().zeros();
            Base::getDerived().grads().zeros();
        }
        else
            check(cudaMemsetAsync(data(), 0, Base::getLength() * sizeof(T)));
    }

#ifdef PHYSICA_HDF5
    template<class Derived>
    auto device_obj<ContinuousVector<Derived>>::read(const H5Loc& loc, const char* name) -> const DataSetType {
        VectorND<T> buffer{};
        auto dataset = buffer.read(loc, name);
        buffer.toDeviceAsync(*this);
        return dataset;
    }

    template<class Derived>
    auto device_obj<ContinuousVector<Derived>>::write(H5Loc& loc, const char* name) const -> DataSetType {
        VectorND<T> buffer{};
        toHost(buffer);
        return buffer.write(loc, name);
    }
#endif

    template<class Derived>
    template<Vector V>
    void ContinuousVector<Derived>::toDevice(device_obj<ContinuousVector<V>>& obj) const {
        toDeviceAsync(obj);
        if constexpr (!std::is_trivially_copy_constructible<V>::value)
            CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector V>
    void ContinuousVector<Derived>::toDeviceAsync(device_obj<ContinuousVector<V>>& obj) const {
        static_assert(std::is_same<T, typename V::ScalarType>::value,
                "[Error]: ScalarType inconsistent, additional buffer is necessary");
        const size_t length = Base::getLength();
        const size_t size = length * sizeof(T);
        obj.resize(length);
        if constexpr (V::SizeAtCompile != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<V>) {
            Base::getDerived().values().toDeviceAsync(obj.getDerived().values());
            Base::getDerived().grads().toDeviceAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyHostToDevice, CUDAContext::getInstance()));
    }
}
