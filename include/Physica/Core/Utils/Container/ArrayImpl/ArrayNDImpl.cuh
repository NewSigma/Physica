/*
 * Copyright 2026 Weibo He.
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

#include "../ArrayND.cuh"

namespace Physica {
    template<class T, int... Dims>
    __host__ __device__ device_obj<ArrayND<T, Dims...>>::device_obj(IndexType shape_, auto&&... args) {
        resize(std::move(shape_), std::forward<decltype(args)>(args)...);
    }

    template<class T, int... Dims>
    __host__ __device__ device_obj<ArrayND<T, Dims...>>::device_obj(std::integral auto... dims) {
        resize(dims...);
    }

    template<class T, int... Dims>
    device_obj<ArrayND<T, Dims...>>::device_obj(const host_obj& storage)
            : arr(storage.arr), shape(storage.shape) {}

    template<class T, int... Dims>
    __device__ auto& device_obj<ArrayND<T, Dims...>>::operator[](this auto&& self, const IndexType& indices) noexcept {
        return *self.data_ptr(indices);
    }

    template<class T, int... Dims>
    __device__ auto& device_obj<ArrayND<T, Dims...>>::operator[](this auto&& self, std::integral auto... dims) noexcept {
        return self[IndexType({static_cast<size_t>(dims)...})];
    }

    template<class T, int... Dims>
    __host__ __device__ void device_obj<ArrayND<T, Dims...>>::resize(IndexType shape_, auto&&... args) {
        if constexpr (IsDevice())
            assert(SizeAtCompile != Dynamic && "[Error]: Do not allocate dynamic tensor in device code");
        else {
            arr.resize(host_obj::toSize(shape_), std::forward<decltype(args)>(args)...);
            if constexpr (StaticShape)
                assert(std::ranges::equal(shape_, IndexType{Dims...}) && "[Error]: Inconsistent size");
            else
                shape = std::move(shape_);
        }
    }

    template<class T, int... Dims>
    __host__ __device__ void device_obj<ArrayND<T, Dims...>>::resize(std::integral auto... dims) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        resize(IndexType{static_cast<size_t>(dims)...});
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::reserve(size_t size) {
        arr.reserve(size);
    }

    template<class T, int... Dims>
    auto device_obj<ArrayND<T, Dims...>>::toHost() const -> host_obj {
        host_obj result;
        toHost(result);
        return result;
    }

    template<class T, int... Dims>
    auto device_obj<ArrayND<T, Dims...>>::toHostAsync() const -> host_obj {
        host_obj result;
        toHostAsync(result);
        return result;
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::toHostAsync(host_obj& obj) const {
        arr.toHostAsync(obj.arr);
        if constexpr (!StaticShape)
            obj.shape = shape;
    }

    template<class T, int... Dims>
    __host__ __device__ size_t device_obj<ArrayND<T, Dims...>>::toIndex1D(const IndexType& indices) const noexcept {
        return IndexType::toIndex1D(getShape(), indices);
    }

    template<class T, int... Dims>
    __host__ __device__ auto device_obj<ArrayND<T, Dims...>>::toIndexND(size_t index) const noexcept -> IndexType {
        assert(index < getSize() && "[Error]: Index out of range");
        return IndexType::toIndexND(getShape(), index);
    }

    template<class T, int... Dims>
    __device__ void device_obj<ArrayND<T, Dims...>>::forND(std::invocable<T&, IndexType> auto func) {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int... Dims>
    __device__ void device_obj<ArrayND<T, Dims...>>::forND(std::invocable<const T&, IndexType> auto func) const {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::zeros() {
        arr.zeros();
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::junk() {
        arr.junk();
    }

    template<class T, int... Dims>
    void device_obj<ArrayND<T, Dims...>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        shape.swap(obj.shape);
    }

    template<class T, int... Dims>
    __host__ __device__ constexpr size_t device_obj<ArrayND<T, Dims...>>::dim(int index) const noexcept {
        assert(0 <= index && index < NDim);
        if constexpr (StaticShape) {
            std::array<int, NDim> buffer{Dims...};
            return buffer[index];
        }
        else
            return shape[index];
    }

    template<class T, int... Dims>
    __host__ __device__ auto device_obj<ArrayND<T, Dims...>>::getShape() const noexcept -> IndexType {
        if constexpr (StaticShape)
            return IndexType({static_cast<size_t>(Dims)...});
        else
            return shape;
    }

    template<class T, int... Dims>
    __host__ __device__ auto* device_obj<ArrayND<T, Dims...>>::data(this auto&& self) noexcept {
        return self.arr.data();
    }

    template<class T, int... Dims>
    __host__ __device__ auto* device_obj<ArrayND<T, Dims...>>::data_ptr(this auto&& self, const IndexType& indices) noexcept {
        return self.data() + self.toIndex1D(indices);
    }

    template<class T, int... Dims>
    __host__ __device__ auto&& device_obj<ArrayND<T, Dims...>>::asArray(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), ArrayType&&>(self.arr);
    }

    template<class T, int... Dims>
    __host__ __device__ size_t device_obj<ArrayND<T, Dims...>>::getSize() const noexcept {
        return arr.getLength();
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::toDeviceAsync() const {
        device_obj<This> result;
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::toDevice(device_obj<This>& obj) const {
        arr.toDevice(obj.arr);
        if constexpr (!StaticShape)
            obj.shape = shape;
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::toDeviceAsync(device_obj<This>& obj) const {
        arr.toDeviceAsync(obj.arr);
        if constexpr (!StaticShape)
            obj.shape = shape;
    }
}
