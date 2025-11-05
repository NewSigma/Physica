/*
 * Copyright 2025 Weibo He.
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

#include "../ArrayND.h"

namespace Physica {
    template<class T, int... Dims>
    ArrayND<T, Dims...>::ArrayND(IndexType shape_, auto&&... args) {
        resize(std::move(shape_), std::forward<decltype(args)>(args)...);
    }

    template<class T, int... Dims>
    ArrayND<T, Dims...>::ArrayND(size_t dim0, auto... dims) {
        resize(dim0, dims...);
    }

    template<class T, int... Dims>
    T& ArrayND<T, Dims...>::operator()(const IndexType& indices) {
        return *data_ptr(indices);
    }

    template<class T, int... Dims>
    const T& ArrayND<T, Dims...>::operator()(const IndexType& indices) const {
        return const_cast<This&>(*this)(indices);
    }

    template<class T, int... Dims>
    T& ArrayND<T, Dims...>::operator()(size_t dim0, auto... dims) {
        return operator()(IndexType({dim0, static_cast<size_t>(dims)...}));
    }

    template<class T, int... Dims>
    const T& ArrayND<T, Dims...>::operator()(size_t dim0, auto... dims) const {
        return const_cast<This&>(*this)(dim0, dims...);
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
    void ArrayND<T, Dims...>::resize(size_t dim0, auto... dims) {
        resize(IndexType({dim0, static_cast<size_t>(dims)...}));
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::zeros() noexcept {
        arr.zeros();
    }

    template<class T, int... Dims>
    size_t ArrayND<T, Dims...>::toIndex1D(const IndexType& indices) const noexcept {
        return toIndex1D(getShape(), indices);
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::toIndexND(size_t index) const noexcept -> IndexType {
        assert(index < getSize() && "[Error]: Index out of range");
        return toIndexND(getShape(), index);
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::forND(std::invocable<T&, IndexType> auto func) {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::forND(std::invocable<const T&, IndexType> auto func) const {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        shape.swap(obj.shape);
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
    auto ArrayND<T, Dims...>::getShape() const noexcept -> IndexType {
        if constexpr (StaticShape)
            return IndexType(Dims...);
        else
            return shape;
    }

    template<class T, int... Dims>
    T* ArrayND<T, Dims...>::data_ptr(const IndexType& indices) noexcept {
        return arr.data() + toIndex1D(indices);
    }

    template<class T, int... Dims>
    const T* ArrayND<T, Dims...>::data_ptr(const IndexType& indices) const noexcept {
        return const_cast<This&>(*this).data_ptr(indices);
    }

    template<class T, int... Dims>
    size_t ArrayND<T, Dims...>::toSize(const IndexType& shape) noexcept {
        const int dim = shape.getLength();
        size_t size = shape[0];
        for (int i = 1; i < dim; ++i)
            size *= shape[i];
        return size;
    }

    template<class T, int... Dims>
    size_t ArrayND<T, Dims...>::toIndex1D(const IndexType& shape, const IndexType& indices) noexcept {
        size_t index = 0;
        size_t stride = 1;
        for (int i = static_cast<int>(shape.getLength()) - 1; i >= 0; --i) {
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            index += indices[i] * stride;
            stride *= shape[i];
        }
        return index;
    }

    template<class T, int... Dims>
    auto ArrayND<T, Dims...>::toIndexND(const IndexType& shape, size_t index) noexcept -> IndexType {
        const int dim = shape.getLength();
        IndexType indices(shape.size());
        size_t remaining = index;
        for (int i = 0; i < dim; ++i) {
            size_t stride = 1;
            for (int j = i + 1; j < dim; ++j)
                stride *= shape[j]; 
            indices[i] = remaining / stride;
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            remaining %= stride;
        }
        return indices;
    }
}
