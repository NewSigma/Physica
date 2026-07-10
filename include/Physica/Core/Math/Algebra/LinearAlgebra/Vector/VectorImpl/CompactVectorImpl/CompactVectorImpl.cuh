/*
 * Copyright 2024-2026 Weibo He.
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

#include "../CompactVector.cuh"

namespace Physica {
    template<class Derived>
    auto device_obj<CompactVector<Derived>>::operator=(Scalar auto x) -> device_obj<Derived>& {
        if constexpr (Base::getSizeAtCompile() == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    void device_obj<CompactVector<Derived>>::reverse(const auto& grad) const noexcept {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename T::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        static_assert(isReverseDiff());
        if constexpr (Scalar<U>)
            Base::getConstCastDerived().grads() += grad;
        else {
            static_assert(DeviceObj<U>, "[Error]: Cannot pass host grad to device");
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
    void device_obj<CompactVector<Derived>>::toHost(CompactVector<V>& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector V>
    void device_obj<CompactVector<Derived>>::toHostAsync(CompactVector<V>& obj) const {
        static_assert(std::same_as<T, typename V::ScalarType>, "[Error]: Incompatible ScalarType");
        const size_t length = Base::getLength();
        const size_t size = length * sizeof(T);
        obj.resize(length);
        if constexpr (obj.getSizeAtCompile() != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<T>) {
            Base::getDerived().values().toHostAsync(obj.getDerived().values());
            Base::getDerived().grads().toHostAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class Derived>
    template<int Size>
    __device__ auto device_obj<CompactVector<Derived>>::packet(size_t index) const noexcept {
        SIMD<T, Size> packet{};
        if constexpr (isReverseDiff())
            packet.load(Base::data_ptr(index).value_ptr());
        else
            packet.load(Base::data_ptr(index));
        return packet;
    }

    template<class Derived>
    template<int Size>
    __device__ auto device_obj<CompactVector<Derived>>::packet(size_t index, size_t count) const noexcept {
        assert(0 < count && count < Size && "[Error]: Invalid size for partial operation");
        assert(index + count <= Base::getLength());
        SIMD<T, Size> packet{};
        if constexpr (isReverseDiff())
            packet.load(Base::data_ptr(index).value_ptr(), count);
        else
            packet.load(Base::data_ptr(index), count);
        return packet;
    }

    template<class Derived>
    __device__ void device_obj<CompactVector<Derived>>::writePacket(const Packet auto packet, size_t index) noexcept {
        using T1 = std::conditional<isReverseDiff(), Tv, T>::type;
        using LocalPacket = std::conditional<packet.size() == 1, T1, SIMD<T1, packet.size()>>::type;
        if constexpr (isReverseDiff())
            LocalPacket(packet).store(Base::data_ptr(index).value_ptr());
        else
            LocalPacket(packet).store(Base::data_ptr(index));
    }

    template<class Derived>
    __device__ void device_obj<CompactVector<Derived>>::writePacket(const Packet auto packet, size_t index, size_t count) noexcept {
        assert(0 < count && count < packet.size() && "[Error]: Invalid size for partial operation");
        assert(index + count <= Base::getLength());
        using T1 = std::conditional<isReverseDiff(), Tv, T>::type;
        using LocalPacket = std::conditional<packet.size() == 1, T1, SIMD<T1, packet.size()>>::type;
        if constexpr (isReverseDiff())
            LocalPacket(packet).store(Base::data_ptr(index).value_ptr(), count);
        else
            LocalPacket(packet).store(Base::data_ptr(index), count);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::head(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<CompactVectorBlock<V, Length>>(std::forward<Self>(self), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::tail(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<CompactVectorBlock<V, Length>>(std::forward<Self>(self), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::segment(this auto&& self, size_t from, size_t to) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<CompactVectorBlock<V, Length>>(std::forward<Self>(self), from, to);
    }

    template<class Derived>
    void device_obj<CompactVector<Derived>>::zeros() {
        if constexpr (Diffable<T>) {
            Base::getDerived().values().zeros();
            Base::getDerived().zero_grad();
        }
        else
            check(cudaMemsetAsync(data(), 0, Base::getLength() * sizeof(T)));
    }

    template<class Derived>
    template<RNG R>
    void device_obj<CompactVector<Derived>>::random_uniform() {
        if constexpr (R::cuRAND_Ready) {
            auto& rng = R::getInstance();
            check(curandSetStream(rng, CUDAContext::getInstance()));
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check(curandGenerateUniform(rng, (float*)data(), length));
            else if constexpr (T::Prec == Float64)
                check(curandGenerateUniformDouble(rng, (double*)data(), length));
            else
                Base::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<CompactVector<Derived>>::random_normal() {
        if constexpr (R::cuRAND_Ready) {
            auto& rng = R::getInstance();
            check(curandSetStream(rng, CUDAContext::getInstance()));
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check(curandGenerateNormal(rng, (float*)data(), length, 0, 1));
            else if constexpr (T::Prec == Float64)
                check(curandGenerateNormalDouble(rng, (double*)data(), length, 0, 1));
            else
                Base::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }

#ifdef PHYSICA_HDF5
    template<class Derived>
    auto device_obj<CompactVector<Derived>>::read(const H5Loc& loc, const char* name) -> const DataSetType {
        VectorND<T> buffer{};
        auto dataset = buffer.read(loc, name);
        buffer.toDeviceAsync(*this);
        return dataset;
    }

    template<class Derived>
    auto device_obj<CompactVector<Derived>>::write(H5Loc& loc, const char* name) const -> DataSetType {
        VectorND<T> buffer{};
        toHost(buffer);
        return buffer.write(loc, name);
    }
#endif

    template<class Derived>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::data_handle(this auto&& self) noexcept {
        return self.data();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactVector<Derived>>::data_ptr(this auto&& self, size_t index) noexcept {
        return self.data_handle() + index;
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<CompactVector<Derived>>::isFastPacket() noexcept {
        return true;
    }

    template<class Derived>
    template<Vector V>
    void CompactVector<Derived>::toDevice(device_obj<CompactVector<V>>& obj) const {
        toDeviceAsync(obj);
        if constexpr (!std::is_trivially_copy_constructible<V>::value)
            CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Vector V>
    void CompactVector<Derived>::toDeviceAsync(device_obj<CompactVector<V>>& obj) const {
        static_assert(std::is_same<T, typename V::ScalarType>::value,
                "[Error]: Type inconsistent between source and target, please cast instead of memcpy");
        const size_t length = Base::getLength();
        const size_t size = length * sizeof(T);
        obj.resize(length);
        if constexpr (obj.getSizeAtCompile() != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<V>) {
            Base::getDerived().values().toDeviceAsync(obj.getDerived().values());
            Base::getDerived().grads().toDeviceAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyHostToDevice, CUDAContext::getInstance()));
    }
}
