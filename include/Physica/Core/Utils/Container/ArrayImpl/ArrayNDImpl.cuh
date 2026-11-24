/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../ArrayND.cuh"

namespace Physica {
    template<class T, int... Dims>
    __host__ __device__ device_obj<ArrayND<T, Dims...>>::device_obj(IndexType shape_, auto&&... args) {
        resize(std::move(shape_), std::forward<decltype(args)>(args)...);
    }

    template<class T, int... Dims>
    __host__ __device__ device_obj<ArrayND<T, Dims...>>::device_obj(std::same_as<size_t> auto... dims) {
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
    __device__ auto& device_obj<ArrayND<T, Dims...>>::operator[](this auto&& self, std::same_as<size_t> auto... dims) noexcept {
        return self[IndexType({dims...})];
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
    __host__ __device__ void device_obj<ArrayND<T, Dims...>>::resize(std::same_as<size_t> auto... dims) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        resize(IndexType{dims...});
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
    __device__ void device_obj<ArrayND<T, Dims...>>::forND(this auto&& self, std::invocable<T&, IndexType> auto func) {
        for (size_t i = 0; i < self.getSize(); ++i)
            func(self.arr[i], self.toIndexND(i));
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
    __host__ __device__ constexpr auto device_obj<ArrayND<T, Dims...>>::getShape() const noexcept -> IndexType {
        if constexpr (StaticShape)
            return IndexType({static_cast<size_t>(Dims)...});
        else
            return shape;
    }

    template<class T, int... Dims>
    __host__ __device__ constexpr size_t device_obj<ArrayND<T, Dims...>>::dim(int index) const noexcept {
        assert(0 <= index && index < NDim);
        if constexpr (StaticShape) {
            Array<int, NDim> buffer{Dims...};
            return buffer[index];
        }
        else
            return shape[index];
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
    __host__ __device__ constexpr size_t device_obj<ArrayND<T, Dims...>>::getSize() const noexcept {
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
