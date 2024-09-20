/*
 * Copyright 2021-2024 Weibo He.
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

#include "DenseMatrixStorage.h"

namespace Physica::Core {
    namespace Internal {
        template<bool IsScalar, class T, size_t Size>
        struct HalfDenseMatrixStorageHelper {
            using Type = Vector<T, Size>;
        };

        template<class T, size_t Size>
        struct HalfDenseMatrixStorageHelper<false, T, Size> {
            using Type = Array<T, Size>;
        };
    }
    /**
     * \class HalfDenseMatrixStorage stores half of the elements of a matrix, while the other half may be symmetric, hermitian, or etc.
     *
     * TODO: This class should extends ArrayBase after merging Core and Utils
     */
    template<class T, size_t Order = Dynamic>
    class HalfDenseMatrixStorage {
        using This = HalfDenseMatrixStorage<T, Order>;
        constexpr static size_t Size = Order * (Order + 1) / 2;
    public:
        using ArrayType = typename Internal::HalfDenseMatrixStorageHelper<is_scalar<T>::value, T, Size>::Type;
        using pointer = typename ArrayType::pointer;
        using const_pointer = typename ArrayType::const_pointer;
        using lvalue_reference = typename ArrayType::lvalue_reference;
        using const_lvalue_reference = typename ArrayType::const_lvalue_reference;
    private:
        ArrayType arr;
        size_t order;
    public:
        HalfDenseMatrixStorage() : arr(), order(Order) {}
        HalfDenseMatrixStorage(size_t order_) : arr(order_ * (order_ + 1) / 2), order(order_) {}
        HalfDenseMatrixStorage(size_t order_, const T& t) : arr(order_ * (order_ + 1) / 2, t), order(order_) {}
        HalfDenseMatrixStorage(size_t row, size_t column);
        HalfDenseMatrixStorage(size_t row, size_t column, const T& t);
        HalfDenseMatrixStorage(std::initializer_list<T> list) : arr(list) {}
        HalfDenseMatrixStorage(const This&) = default;
        HalfDenseMatrixStorage(This&&) noexcept = default;
        ~HalfDenseMatrixStorage() = default;
        /* Operators */
        HalfDenseMatrixStorage& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ lvalue_reference operator[](size_t index) { return arr[index]; }
        [[nodiscard]] __host__ __device__ const_lvalue_reference operator[](size_t index) const { return arr[index]; }
        [[nodiscard]] __host__ __device__ inline lvalue_reference operator()(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline const_lvalue_reference operator()(size_t row, size_t col) const;
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t column, Args&&... args);
        void resize(size_t row) { resize(row, row); }
        void swap(This& __restrict storage) noexcept;
        /* Getters */
        [[nodiscard]] const ArrayType& getArray() const noexcept { return arr; }
        [[nodiscard]] ArrayType& getArray() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return order; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return arr.size(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return arr.getCapacity(); }
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return arr.data(); }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return arr.data(); }
        [[nodiscard]] __host__ __device__ inline T* data_ptr(size_t row, size_t column) noexcept;
        [[nodiscard]] __host__ __device__ inline const T* data_ptr(size_t row, size_t column) const noexcept;
        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const noexcept;
    };

    template<class T, size_t Order>
    HalfDenseMatrixStorage<T, Order>::HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column)
            : HalfDenseMatrixStorage(row) {
        assert(row == column);
    }

    template<class T, size_t Order>
    HalfDenseMatrixStorage<T, Order>::HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column, const T& t)
            : HalfDenseMatrixStorage(row, t) {
        assert(row == column);
    }

    template<class T, size_t Order>
    __host__ __device__ inline typename HalfDenseMatrixStorage<T, Order>::lvalue_reference
    HalfDenseMatrixStorage<T, Order>::operator()(size_t row, size_t col) { return (*this)[toIndex1D(row, col)]; }

    template<class T, size_t Order>
    __host__ __device__ inline typename HalfDenseMatrixStorage<T, Order>::const_lvalue_reference
    HalfDenseMatrixStorage<T, Order>::operator()(size_t row, size_t col) const { return (*this)[toIndex1D(row, col)]; }

    template<class T, size_t Order>
    template<class... Args>
    void HalfDenseMatrixStorage<T, Order>::resize(size_t row, [[maybe_unused]] size_t column, Args&&... args) {
        assert(row == column);
        arr.resize(row * (row + 1) / 2, std::forward<Args>(args)...);
        order = row;
    }

    template<class T, size_t Order>
    void HalfDenseMatrixStorage<T, Order>::swap(This& __restrict storage) noexcept {
        assert(this != &storage && "[Error]: Self swap is likely a bug");
        arr.swap(storage.arr);
        std::swap(order, storage.order);
    }

    template<class T, size_t Order>
    __host__ __device__ inline T* HalfDenseMatrixStorage<T, Order>::data_ptr(size_t row, size_t column) noexcept {
        return arr.data() + toIndex1D(row, column);
    }

    template<class T, size_t Order>
    __host__ __device__ inline const T* HalfDenseMatrixStorage<T, Order>::data_ptr(size_t row, size_t column) const noexcept {
        return arr.data() + toIndex1D(row, column);
    }

    template<class T, size_t Order>
    __host__ __device__ size_t HalfDenseMatrixStorage<T, Order>::toIndex1D(size_t r, size_t c) const noexcept {
        const size_t order = getOrder();
        assert(r < order && c < order);
        const bool exchange = c < r;
        const size_t min = exchange ? c : r;
        const size_t max = exchange ? r : c;
        return (order * 2U - min - 1) * min / 2U + max;
    }
}

/*namespace Physica {
    template<class T, size_t Order>
    class Traits<Core::HalfDenseMatrixStorage<T, Order>> {
        constexpr static bool IsScalar = Core::is_scalar<T>::value;
        constexpr static size_t Size = Order * (Order + 1) / 2;
    public:
        using ArrayType = typename std::conditional<IsScalar, Core::Vector<T, Size>, Core::Array<T, Size>>::type;
    };
}*/

namespace std {
    template<class T, size_t Order>
    inline void swap(Physica::Core::HalfDenseMatrixStorage<T, Order>& __restrict mat1,
                     Physica::Core::HalfDenseMatrixStorage<T, Order>& __restrict mat2) noexcept {
        mat1.swap(mat2);
    }
}
