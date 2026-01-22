/*
 * Copyright 2025-2026 Weibo He.
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

#include "DiagMatrix.h"
#include "MatrixImpl/RValueMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"

namespace Physica {
    template<Scalar T, size_t Order>
    class device_obj<DiagMatrix<T, Order>> : public device_obj<RValueMatrix<DiagMatrix<T, Order>>> {
        using host_obj = DiagMatrix<T, Order>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using VectorType = device_obj<DenseVector<T, Order>>;
    protected:
        using typename Base::Tr;
        using typename Base::Trv;
        using typename Base::Tv;
    private:
        VectorType diags;
    public:
        device_obj() = default;
        explicit device_obj(size_t order);
        explicit device_obj(const Vector auto& diags_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        __host__ __device__ void operator*=(const Scalar auto& x);
        [[nodiscard]] __host__ __device__ auto operator*(Vector auto&& v) const noexcept;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] This inv() const;
        [[nodiscard]] __host__ __device__ const This& transpose() const noexcept { return *this; }

        void resize(size_t order);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& diag(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getOrder(); }
        /* Static members */
        [[nodiscard]] static This identity(size_t order);
    };

    template<Scalar T, size_t Order>
    device_obj<DiagMatrix<T, Order>>::device_obj(size_t order) : diags(order) {}

    template<Scalar T, size_t Order>
    device_obj<DiagMatrix<T, Order>>::device_obj(const Vector auto& diags_) : diags(diags_) {}

    template<Scalar T, size_t Order>
    __host__ __device__ void device_obj<DiagMatrix<T, Order>>::operator*=(const Scalar auto& x) {
        diags *= x;
    }

    template<Scalar T, size_t Order>
    __host__ __device__ auto device_obj<DiagMatrix<T, Order>>::operator*(Vector auto&& v) const noexcept {
        assert(getCol() == v.getLength());
        return hadamard(diags, std::forward<decltype(v)>(v));
    }

    template<Scalar T, size_t Order>
    __device__ T device_obj<DiagMatrix<T, Order>>::calc(size_t row, size_t col) const {
        return row == col ? diags[row] : T(0);
    }

    template<Scalar T, size_t Order>
    __device__ auto device_obj<DiagMatrix<T, Order>>::calc_value(size_t row, size_t col) const -> Tv {
        return calc(row, col).value();
    }

    template<Scalar T, size_t Order>
    auto device_obj<DiagMatrix<T, Order>>::inv() const -> This {
        return This(reciprocal(diags));
    }

    template<Scalar T, size_t Order>
    void device_obj<DiagMatrix<T, Order>>::resize(size_t order) {
        diags.resize(order);
    }

    template<Scalar T, size_t Order>
    void device_obj<DiagMatrix<T, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        diags.swap(obj.diags);
    }

    template<Scalar T, size_t Order>
    __host__ __device__ auto&& device_obj<DiagMatrix<T, Order>>::diag(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.diags);
    }

    template<Scalar T, size_t Order>
    __host__ __device__ size_t device_obj<DiagMatrix<T, Order>>::getOrder() const noexcept {
        if constexpr (Order == Dynamic)
            return diags.getLength();
        return Order;
    }

    template<Scalar T, size_t Order>
    auto device_obj<DiagMatrix<T, Order>>::identity(size_t order) -> This {
        return This(VectorType(order, 1));
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<device_obj<DiagMatrix<T, Order>>> : public Traits<DiagMatrix<T, Order>> {};
}
