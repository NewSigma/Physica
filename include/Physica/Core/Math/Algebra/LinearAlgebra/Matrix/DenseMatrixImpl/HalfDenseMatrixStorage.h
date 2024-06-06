/*
 * Copyright 2021-2022 WeiBo He.
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

namespace Physica::Core::Internal {
    template<class T, size_t Order, size_t MaxOrder> class HalfDenseMatrixStorage;

    template<class T, size_t Order, size_t MaxOrder>
    class Traits<HalfDenseMatrixStorage<T, Order, MaxOrder>> {
        constexpr static bool IsScalar = is_scalar<T>::value;
        constexpr static size_t Size = Order * (Order + 1) / 2;
        constexpr static size_t MaxSize = MaxOrder * (MaxOrder + 1) / 2;
    public:
        using Base = typename std::conditional<IsScalar, Vector<T, Size, MaxSize>, Utils::Array<T, Size, MaxSize>>::type;
    };

    template<class T, size_t Order, size_t MaxOrder>
    class HalfDenseMatrixStorage : public Traits<HalfDenseMatrixStorage<T, Order, MaxOrder>>::Base {
        using This = HalfDenseMatrixStorage<T, Order, MaxOrder>;
    public:
        using Base = typename Traits<This>::Base;
    private:
        size_t order;
    public:
        HalfDenseMatrixStorage() : Base(), order(MaxOrder) {}
        HalfDenseMatrixStorage(size_t order_) : Base(order_ * (order_ + 1) / 2), order(order_) {}
        HalfDenseMatrixStorage(size_t order_, const T& t) : Base(order_ * (order_ + 1) / 2, t), order(order_) {}
        HalfDenseMatrixStorage(size_t row, size_t column);
        HalfDenseMatrixStorage(size_t row, size_t column, const T& t);
        HalfDenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        HalfDenseMatrixStorage(const This&) = default;
        HalfDenseMatrixStorage(This&&) noexcept = default;
        ~HalfDenseMatrixStorage() = default;
        /* Operators */
        HalfDenseMatrixStorage& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void resize(size_t order_) { Base::resize(order_ * (order_ + 1) / 2); order = order_; }
        void resize(size_t row, [[maybe_unused]] size_t column) { assert(row == column); resize(row); order = row; } //Necessary to CRTP
        void swap(HalfDenseMatrixStorage& __restrict storage) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return order; }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ inline T* data_ptr(size_t row, size_t column) noexcept;
        [[nodiscard]] __host__ __device__ inline const T* data_ptr(size_t row, size_t column) const noexcept;
        [[nodiscard]] size_t toIndex1D(size_t r, size_t c) const noexcept;
    };

    template<class T, size_t Order, size_t MaxOrder>
    HalfDenseMatrixStorage<T, Order, MaxOrder>::HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column)
            : HalfDenseMatrixStorage(row) {
        assert(row == column);
    }

    template<class T, size_t Order, size_t MaxOrder>
    HalfDenseMatrixStorage<T, Order, MaxOrder>::HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column, const T& t)
            : HalfDenseMatrixStorage(row, t) {
        assert(row == column);
    }

    template<class T, size_t Order, size_t MaxOrder>
    void HalfDenseMatrixStorage<T, Order, MaxOrder>::swap(HalfDenseMatrixStorage<T, Order, MaxOrder>& __restrict storage) noexcept {
        assert(this != &storage && "[Error]: Self swap is likely a bug");
        std::swap(order, storage.order);
        Base::swap(storage);
    }

    template<class T, size_t Order, size_t MaxOrder>
    __host__ __device__ inline T* HalfDenseMatrixStorage<T, Order, MaxOrder>::data_ptr(size_t row, size_t column) noexcept {
        return Base::data() + toIndex1D(row, column);
    }

    template<class T, size_t Order, size_t MaxOrder>
    __host__ __device__ inline const T* HalfDenseMatrixStorage<T, Order, MaxOrder>::data_ptr(size_t row, size_t column) const noexcept {
        return Base::data() + toIndex1D(row, column);
    }

    template<class T, size_t Order, size_t MaxOrder>
    size_t HalfDenseMatrixStorage<T, Order, MaxOrder>::toIndex1D(size_t r, size_t c) const noexcept {
        const size_t order = getOrder();
        assert(r < order && c < order);
        const bool exchange = c < r;
        const size_t min = exchange ? c : r;
        const size_t max = exchange ? r : c;
        return (order * 2U - min - 1) * min / 2U + max;
    }

    template<class T, size_t Order, size_t MaxOrder>
    inline void swap(HalfDenseMatrixStorage<T, Order, MaxOrder>& __restrict mat1,
                     HalfDenseMatrixStorage<T, Order, MaxOrder>& __restrict mat2) noexcept {
        mat1.swap(mat2);
    }
}
