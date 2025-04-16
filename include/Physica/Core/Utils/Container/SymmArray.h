/*
 * Copyright 2021-2025 Weibo He.
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
    class SymmArray : public ArrayBase<SymmArray<T, Order>, HostAllocator<T>> {
        using This = SymmArray<T, Order>;
        using Base = ArrayBase<This, HostAllocator<T>>;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
    protected:
        using ArrayType = Traits<This>::ArrayType;
    private:
        ArrayType arr;
        size_t order;
    public:
        SymmArray() : arr(), order(Order) {}
        SymmArray(size_t order_);
        SymmArray(size_t order_, const T& t);
        SymmArray(size_t row, size_t col);
        SymmArray(size_t row, size_t col, const T& t);
        SymmArray(std::initializer_list<T> list) : arr(list) {}
        SymmArray(const This&) = default;
        SymmArray(This&&) noexcept = default;
        ~SymmArray() = default;
        /* Operators */
        SymmArray& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ inline lvalue_reference operator()(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline const_lvalue_reference operator()(size_t row, size_t col) const;
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t col, Args&&... args);
        void resize(size_t row) { resize(row, row); }
        void swap(This& __restrict storage) noexcept;

        inline const H5DataSet<1> read(const H5Loc& loc, const char* name);
        inline H5DataSet<1> write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return order; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return arr.size(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return arr.getCapacity(); }
        [[nodiscard]] __host__ __device__ auto data() noexcept { return arr.data(); }
        [[nodiscard]] __host__ __device__ auto data() const noexcept { return arr.data(); }
        [[nodiscard]] __host__ __device__ inline auto data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline auto data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const noexcept;
    };

    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t order_) : arr(order_ * (order_ + 1) / 2), order(order_) {}
    /**
     * If std::same_as<T, size_t> is true, the semantics is ambiguous.
     *
     * Use either SymmArray(size_t) or SymmArray(size_t, size_t, size_t) in this case.
     */
    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t order_, const T& t)
            : arr(order_ * (order_ + 1) / 2, t), order(order_) {}

    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t row, [[maybe_unused]] size_t col)
            : SymmArray(row) {
        assert(row == col);
    }

    template<class T, size_t Order>
    SymmArray<T, Order>::SymmArray(size_t row, [[maybe_unused]] size_t col, const T& t)
            : SymmArray(row, t) {
        assert(row == col);
    }

    template<class T, size_t Order>
    __host__ __device__ inline SymmArray<T, Order>::lvalue_reference
    SymmArray<T, Order>::operator()(size_t row, size_t col) {
        return (*this)[toIndex1D(row, col)];
    }

    template<class T, size_t Order>
    __host__ __device__ inline SymmArray<T, Order>::const_lvalue_reference
    SymmArray<T, Order>::operator()(size_t row, size_t col) const {
        return (*this)[toIndex1D(row, col)];
    }

    template<class T, size_t Order>
    template<class... Args>
    void SymmArray<T, Order>::resize(size_t row, [[maybe_unused]] size_t col, Args&&... args) {
        assert(row == col);
        arr.resize(row * (row + 1) / 2, std::forward<Args>(args)...);
        order = row;
    }

    template<class T, size_t Order>
    void SymmArray<T, Order>::swap(This& __restrict storage) noexcept {
        assert(this != &storage && "[Error]: Self swap is likely a bug");
        arr.swap(storage.arr);
        std::swap(order, storage.order);
    }

#ifdef PHYSICA_HDF5
    template<class T, size_t Order>
    inline const H5DataSet<1> SymmArray<T, Order>::read(const H5Loc& loc, const char* name) {
        auto group = arr.read(loc, name);
        assert(order * (order + 1) / 2 == size() && "[Error]: Order is not well initialized");
        return group;
    }

    template<class T, size_t Order>
    inline H5DataSet<1> SymmArray<T, Order>::write(H5Loc& loc, const char* name) const {
        return arr.write(loc, name);
    }
#endif

    template<class T, size_t Order>
    __host__ __device__ inline auto SymmArray<T, Order>::data_ptr(size_t row, size_t col) noexcept {
        return arr.data() + toIndex1D(row, col);
    }

    template<class T, size_t Order>
    __host__ __device__ inline auto SymmArray<T, Order>::data_ptr(size_t row, size_t col) const noexcept {
        return arr.data() + toIndex1D(row, col);
    }

    template<class T, size_t Order>
    __host__ __device__ size_t SymmArray<T, Order>::toIndex1D(size_t r, size_t c) const noexcept {
        assert(r < order && c < order);
        const bool exchange = c < r;
        const size_t min = exchange ? c : r;
        const size_t max = exchange ? r : c;
        return (order * 2U - min - 1) * min / 2U + max;
    }
}

namespace Physica {
    template<class T, size_t Order>
    class Traits<SymmArray<T, Order>> {
        template<bool, size_t Size>
        struct Helper {
            using Type = DenseVector<T, Size>;
        };

        template<size_t Size>
        struct Helper<false, Size> {
            using Type = Array<T, Size>;
        };

        constexpr static bool IsScalar = Scalar<T>;
        constexpr static size_t SizeAtCompile = Order * (Order + 1) / 2;
    public:
        using ElemType = T;
        using ArrayType = Helper<IsScalar, SizeAtCompile>::Type;
    };
}

namespace std {
    template<class T, size_t Order>
    inline void swap(Physica::SymmArray<T, Order>& __restrict mat1,
                     Physica::SymmArray<T, Order>& __restrict mat2) noexcept {
        mat1.swap(mat2);
    }
}
