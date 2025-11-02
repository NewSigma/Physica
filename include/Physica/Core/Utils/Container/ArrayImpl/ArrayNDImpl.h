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
    template<class T, int Dim>
    ArrayND<T, Dim>::ArrayND(ShapeType shape_, auto&&... args) {
        resize(std::move(shape_), std::forward<decltype(args)>(args)...);
    }

    template<class T, int Dim>
    ArrayND<T, Dim>::ArrayND(size_t dim0, auto... dims) {
        resize(dim0, dims...);
    }

    template<class T, int Dim>
    T& ArrayND<T, Dim>::operator()(const ShapeType& indices) {
        return *data_ptr(indices);
    }

    template<class T, int Dim>
    const T& ArrayND<T, Dim>::operator()(const ShapeType& indices) const {
        return const_cast<This&>(*this)(indices);
    }

    template<class T, int Dim>
    T& ArrayND<T, Dim>::operator()(size_t dim0, auto... dims) {
        return operator()(toShape(dim0, dims...));
    }

    template<class T, int Dim>
    const T& ArrayND<T, Dim>::operator()(size_t dim0, auto... dims) const {
        return const_cast<This&>(*this)(dim0, dims...);
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::resize(ShapeType shape_, auto&&... args) {
        assert(shape_.getLength() < std::numeric_limits<int>::max() && "[Error]: Does not support");
        shape = std::move(shape_);
        arr.resize(toSize(shape), std::forward<decltype(args)>(args)...);
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::resize(size_t dim0, auto... dims) {
        resize(toShape(dim0, dims...));
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::zeros() noexcept {
        arr.zeros();
    }

    template<class T, int Dim>
    size_t ArrayND<T, Dim>::toIndex1D(const ShapeType& indices) const noexcept {
        return toIndex1D(shape, indices);
    }

    template<class T, int Dim>
    auto ArrayND<T, Dim>::toIndexND(size_t index) const noexcept -> ShapeType {
        assert(index < getSize() && "[Error]: Index out of range");
        return toIndexND(shape, index);
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::forND(std::invocable<T&, ShapeType> auto func) {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::forND(std::invocable<const T&, ShapeType> auto func) const {
        for (size_t i = 0; i < getSize(); ++i)
            func(arr[i], toIndexND(i));
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        shape.swap(obj.shape);
    }

    template<class T, int Dim>
    T* ArrayND<T, Dim>::data_ptr(const ShapeType& indices) noexcept {
        return arr.data() + toIndex1D(indices);
    }

    template<class T, int Dim>
    const T* ArrayND<T, Dim>::data_ptr(const ShapeType& indices) const noexcept {
        return const_cast<This&>(*this).data_ptr(indices);
    }

    template<class T, int Dim>
    size_t ArrayND<T, Dim>::toSize(const ShapeType& shape) noexcept {
        const int dim = shape.getLength();
        size_t size = shape[0];
        for (int i = 1; i < dim; ++i)
            size *= shape[i];
        return size;
    }

    template<class T, int Dim>
    size_t ArrayND<T, Dim>::toIndex1D(const ShapeType& shape, const ShapeType& indices) noexcept {
        size_t index = 0;
        size_t stride = 1;
        for (int i = static_cast<int>(shape.getLength()) - 1; i >= 0; --i) {
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            index += indices[i] * stride;
            stride *= shape[i];
        }
        return index;
    }

    template<class T, int Dim>
    auto ArrayND<T, Dim>::toIndexND(const ShapeType& shape, size_t index) noexcept -> ShapeType {
        const int dim = shape.getLength();
        ShapeType indices(shape.size());
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

    template<class T, int Dim>
    auto ArrayND<T, Dim>::toShape(auto... dims) noexcept -> ShapeType {
        constexpr int Dim1 = sizeof...(dims);
        static_assert(Dim == Dim1 || Dim == Dynamic, "[Error]: Dim is not self consistent");
        ShapeType indices(Dim1);
        toShapeImpl(indices, 0, dims...);
        return indices;
    }

    template<class T, int Dim>
    void ArrayND<T, Dim>::toShapeImpl(ShapeType& arr, int count, size_t dim0, auto... dims) noexcept {
        arr[count] = dim0;
        toShapeImpl(arr, count + 1, dims...);
    }
}
