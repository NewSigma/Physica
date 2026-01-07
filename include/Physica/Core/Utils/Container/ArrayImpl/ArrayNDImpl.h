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
    ArrayND<T, Dims...>::ArrayND(std::integral auto... dims) {
        resize(dims...);
    }

    template<class T, int... Dims>
    T& ArrayND<T, Dims...>::operator[](const IndexType& indices) {
        return *data_ptr(indices);
    }

    template<class T, int... Dims>
    const T& ArrayND<T, Dims...>::operator[](const IndexType& indices) const {
        return const_cast<This&>(*this)[indices];
    }

    template<class T, int... Dims>
    T& ArrayND<T, Dims...>::operator[](std::integral auto... dims) {
        return operator[](IndexType({static_cast<size_t>(dims)...}));
    }

    template<class T, int... Dims>
    const T& ArrayND<T, Dims...>::operator[](std::integral auto... dims) const {
        return const_cast<This&>(*this)[dims...];
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
    void ArrayND<T, Dims...>::resize(std::integral auto... dims) {
        resize(IndexType({static_cast<size_t>(dims)...}));
    }

    template<class T, int... Dims>
    void ArrayND<T, Dims...>::zeros() noexcept {
        arr.zeros();
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
    auto* ArrayND<T, Dims...>::data(this auto&& self) noexcept {
        return self.arr.data();
    }

    template<class T, int... Dims>
    auto* ArrayND<T, Dims...>::data_ptr(this auto&& self, const IndexType& indices) noexcept {
        return self.data() + self.toIndex1D(indices);
    }

    template<class T, int... Dims>
    size_t ArrayND<T, Dims...>::toSize(const IndexType& shape) noexcept {
        const int dim = shape.getLength();
        size_t size = shape[0];
        for (int i = 1; i < dim; ++i)
            size *= shape[i];
        return size;
    }
}
