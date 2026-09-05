/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    /**
     * \class SymmArray stores half of the elements of a matrix, while the other half may be symmetric, hermitian, or etc.
     */
    template<class T, size_t Order = Dynamic>
    class SymmArray : public ArrayMixin<SymmArray<T, Order>, HostAllocator<T>> {
        using This = SymmArray<T, Order>;
        using Base = ArrayMixin<This, HostAllocator<T>>;
    protected:
        using ArrayType = Traits<This>::ArrayType;
    private:
        ArrayType arr;
        size_t order = Order;
    public:
        SymmArray() = default;
        SymmArray(size_t order);
        SymmArray(size_t order, const T& t);
        SymmArray(const This&) = default;
        SymmArray(This&&) noexcept = default;
        ~SymmArray() = default;
        /* Operators */
        SymmArray& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ decltype(auto) operator[](this auto&&, size_t row, size_t col);
        /* Operations */
        void resize(size_t order_, auto&&... args);
        void zeros() noexcept { arr.zeros(); }
        void swap(This& __restrict storage) noexcept;

        const H5Dataset<1> read(const H5Loc& loc, const char* name);
        H5Dataset<1> write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] auto&& asArray(this auto&& self) noexcept { return self.arr; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return order; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return arr.getCapacity(); }
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const noexcept;
    };

    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t order) : arr(order * (order + 1) / 2), order(order) {}
    /**
     * If std::same_as<T, size_t> is true, the semantics is ambiguous.
     *
     * Use either SymmArray(size_t) or SymmArray(size_t, size_t, size_t) in this case.
     */
    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t order, const T& t)
            : arr(order * (order + 1) / 2, t), order(order) {}

    template<class T, size_t Order>
    __host__ __device__ decltype(auto) SymmArray<T, Order>::operator[](this auto&& self, size_t row, size_t col) {
        return self.arr[self.toIndex1D(row, col)];
    }

    template<class T, size_t Order>
    void SymmArray<T, Order>::resize(size_t order_, auto&&... args) {
        arr.resize(order_ * (order_ + 1) / 2, std::forward<decltype(args)>(args)...);
        order = order_;
    }

    template<class T, size_t Order>
    void SymmArray<T, Order>::swap(This& __restrict storage) noexcept {
        assert(this != &storage && "[Error]: Self swap is likely a bug");
        arr.swap(storage.arr);
        std::swap(order, storage.order);
    }

#ifdef PHYSICA_HDF5
    template<class T, size_t Order>
    const H5Dataset<1> SymmArray<T, Order>::read(const H5Loc& loc, const char* name) {
        auto group = arr.read(loc, name);
        assert(order * (order + 1) / 2 == getLength() && "[Error]: Order is not well initialized");
        return group;
    }

    template<class T, size_t Order>
    H5Dataset<1> SymmArray<T, Order>::write(H5Loc& loc, const char* name) const {
        return arr.write(loc, name);
    }
#endif

    template<class T, size_t Order>
    __host__ __device__ auto SymmArray<T, Order>::data(this auto&& self) noexcept {
        return self.arr.data();
    }

    template<class T, size_t Order>
    __host__ __device__ auto SymmArray<T, Order>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        return self.data() + self.toIndex1D(row, col);
    }

    template<class T, size_t Order>
    __host__ __device__ size_t SymmArray<T, Order>::toIndex1D(size_t r, size_t c) const noexcept {
        assert(r < order && c < order);
        const bool exchange = c < r;
        const size_t min = exchange ? c : r;
        const size_t max = exchange ? r : c;
        return (order * 2U - min - 1) * min / 2U + max;
    }

    template<class T, size_t Order>
    void swap(SymmArray<T, Order>& mat1, SymmArray<T, Order>& mat2) noexcept {
        mat1.swap(mat2);
    }
}

namespace Physica {
    template<class T, size_t Order>
    class Traits<SymmArray<T, Order>> {
        template<bool, size_t Length>
        struct Helper {
            using Type = DenseVector<T, Length>;
        };

        template<size_t Length>
        struct Helper<false, Length> {
            using Type = Array<T, Length>;
        };
    public:
        using value_type = T;
        using ArrayType = Helper<Scalar<T>, Order * (Order + 1) / 2>::Type;
    };
}
