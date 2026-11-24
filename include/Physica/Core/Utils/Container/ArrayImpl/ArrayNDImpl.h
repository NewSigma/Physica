/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../ArrayND.h"

namespace Physica {
    template<class T, int... Dims>
    ArrayND<T, Dims...>::ArrayND(IndexType shape_, auto&&... args) {
        resize(std::move(shape_), std::forward<decltype(args)>(args)...);
    }

    template<class T, int... Dims>
    ArrayND<T, Dims...>::ArrayND(std::same_as<size_t> auto... dims) {
        resize(dims...);
    }

    template<class T, int... Dims>
    auto& ArrayND<T, Dims...>::operator[](this auto&& self, const IndexType& indices) noexcept {
        return *self.data_ptr(indices);
    }

    template<class T, int... Dims>
    auto& ArrayND<T, Dims...>::operator[](this auto&& self, std::same_as<size_t> auto... dims) noexcept {
        return self[IndexType({dims...})];
    }

    template<class T, int... Dims>
    bool ArrayND<T, Dims...>::operator==(this const This& lhs, const This& rhs) {
        if constexpr (!StaticShape)
            if (lhs.shape != rhs.shape)
                return false;
        return lhs.arr == rhs.arr;
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::resize(IndexType shape_, auto&&... args) {
        arr.resize(toSize(shape_), std::forward<decltype(args)>(args)...);
        if constexpr (StaticShape)
            assert(std::ranges::equal(shape_, IndexType{Dims...}) && "[Error]: Inconsistent size");
        else
            shape = std::move(shape_);
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::resize(std::same_as<size_t> auto... dims) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        resize(IndexType{dims...});
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::reserve(size_t size) noexcept {
        arr.reserve(size);
    }

    template<class T, int... Dims>
    size_t ArrayND<T, Dims...>::toIndex1D(const IndexType& indices) const noexcept {
        return IndexType::toIndex1D(getShape(), indices);
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::toIndexND(size_t index) const noexcept -> IndexType {
        assert(index < getSize() && "[Error]: Index out of range");
        return IndexType::toIndexND(getShape(), index);
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::forND(this auto&& self, std::invocable<T&, IndexType> auto func) {
        Physica::forND(self.getShape(), [&self, func](const IndexType& index) mutable {
            func(*self.data_ptr(index), index);
        });
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::zeros() noexcept {
        arr.zeros();
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::junk() noexcept {
        arr.junk();
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        shape.swap(obj.shape);
    }

    template<class T, int... Dims>
    constexpr auto ArrayND<T, Dims...>::getShape() const noexcept -> IndexType {
        if constexpr (StaticShape)
            return IndexType({static_cast<size_t>(Dims)...});
        else
            return shape;
    }

    template<class T, int... Dims>
    constexpr size_t ArrayND<T, Dims...>::dim(int index) const noexcept {
        assert(0 <= index && index < NDim);
        if constexpr (StaticShape) {
            std::array<int, NDim> buffer{Dims...};
            return buffer[index];
        }
        else
            return shape[index];
    }

    template<class T, int... Dims>
    auto* ArrayND<T, Dims...>::data(this auto&& self) noexcept {
        return self.arr.data();
    }

    template<class T, int... Dims>
    auto* ArrayND<T, Dims...>::data_ptr(this auto&& self, const IndexType& indices) noexcept {
        return self.data() + self.toIndex1D(indices);
    }

    template<class T, int... Dims>
    auto&& ArrayND<T, Dims...>::asArray(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), ArrayType&&>(self.arr);
    }

    template<class T, int... Dims>
    __host__ __device__ size_t ArrayND<T, Dims...>::toSize(const IndexType& shape) noexcept {
        return Index<NDim>::toSize(shape);
    }

    template<class T, int... Dims>
    consteval auto ArrayND<T, Dims...>::getShapeAtCompile() noexcept -> IndexType {
        if constexpr (StaticShape)
            return IndexType({static_cast<size_t>(Dims)...});
        else
            return IndexType(NDim, Dynamic);
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::zeros(IndexType shape_) -> This {
        auto result = This(std::move(shape_));
        result.zeros();
        return result;
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::zeros(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        return zeros(IndexType({dims...}));
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::junk(IndexType shape_) -> This {
        auto result = This(std::move(shape_));
        result.junk();
        return result;
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::junk(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        return junk(IndexType({dims...}));
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::read(IndexType shape_, const T* __restrict p) noexcept -> This {
        This result;
        if constexpr (!StaticShape)
            result.shape = std::move(shape_);
        result.arr = ArrayType::read(toSize(result.getShape()), p);
        return result;
    }
}
