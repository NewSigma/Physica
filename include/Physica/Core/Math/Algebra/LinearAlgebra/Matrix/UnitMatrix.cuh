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


#include "MatrixImpl/RValueMatrix.cuh"
#include "UnitMatrix.h"

namespace Physica {
    template<Scalar T, size_t Order>
    class device_obj<UnitMatrix<T, Order>> : public device_obj<RValueMatrix<UnitMatrix<T, Order>>> {
        using host_obj = UnitMatrix<T, Order>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using IndexType = host_obj::IndexType;
    protected:
        using typename Base::Tv;
    private:
        [[no_unique_address]] IndexType order;
    public:
        device_obj() = default;
        explicit device_obj(size_t order_);
        explicit device_obj(const Matrix auto& m);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        template<Vector V>
        [[nodiscard]] __host__ __device__ V&& operator*(V&& v) const noexcept;
        /* Operations */
        template<Matrix M>
        void assign(M&& target) const;

        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ T calc_value(size_t row, size_t col) const;

        [[nodiscard]] __host__ __device__ const This& transpose() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const This& conjugate() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const This& hermite() const noexcept { return *this; }
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
    };

    template<Scalar T, size_t Order>
    device_obj<UnitMatrix<T, Order>>::device_obj(size_t order_) : order(order_) {
        assert((Order == Dynamic || Order == order_) && "[Error]: tparam and param is not consistent");
    }

    template<Scalar T, size_t Order>
    device_obj<UnitMatrix<T, Order>>::device_obj(const Matrix auto& m) : device_obj(m.getRow()) {
        assert(m.isSquare());
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    void device_obj<UnitMatrix<T, Order>>::assign(M&& target) const {
        target.zeros();
        auto setones = [target_ = asStruct(target)] __device__ () mutable {
            auto& target = target_.getDerived();
            auto i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i < target.getRow())
                target(i, i) = T(1);
        };

        constexpr int WarpSize = CUDADevAttr::WarpSize;
        const int numBlock = (getRow() + (WarpSize - 1)) / WarpSize;
        CUDAExecutor::launch(setones, {numBlock, WarpSize});
    }

    template<Scalar T, size_t Order>
    template<Vector V>
    __host__ __device__ V&& device_obj<UnitMatrix<T, Order>>::operator*(V&& v) const noexcept {
        static_assert(CUDA<V>);
        assert(getCol() == v.getLength());
        return std::forward<V>(v);
    }

    template<Scalar T, size_t Order>
    __device__ T device_obj<UnitMatrix<T, Order>>::calc(size_t row, size_t col) const {
        return T(row == col ? 1 : 0);
    }

    template<Scalar T, size_t Order>
    __device__ T device_obj<UnitMatrix<T, Order>>::calc_value(size_t row, size_t col) const {
        return calc(row, col);
    }

    template<Scalar T, size_t Order>
    __host__ __device__ void device_obj<UnitMatrix<T, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(order, obj.order);
    }

    template<Scalar T, size_t Order>
    __host__ __device__ size_t device_obj<UnitMatrix<T, Order>>::getRow() const noexcept {
        if constexpr (Order == Dynamic)
            return order;
        return Order;
    }

    template<Scalar T, size_t Order>
    __host__ __device__ size_t device_obj<UnitMatrix<T, Order>>::getCol() const noexcept {
        return getRow();
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<device_obj<UnitMatrix<T, Order>>> : public Traits<UnitMatrix<T, Order>> {};
}

#include "UnitMatrixImpl/Mul.h"
