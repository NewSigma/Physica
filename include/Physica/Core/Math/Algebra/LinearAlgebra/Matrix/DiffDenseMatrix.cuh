/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/CompactMatrix.cuh"
#include "DenseMatrix.cuh"
#include "DiffDenseMatrix.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int Major>
    class device_obj<DenseMatrix<Diff<T, Mode, Order>, Major>>
            : public device_obj<CompactMatrix<DenseMatrix<Diff<T, Mode, Order>, Major>>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<device_obj<DenseMatrix<Diff<T, Mode, Order>, Major>>>, Empty>::type {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
        using host_obj = DenseMatrix<Diff<T, Mode, Order>, Major>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactMatrix<host_obj>>;
    public:
        using ScalarType = Base::ScalarType;
        template<Scalar U>
        using rebind_scalar = device_obj<DenseMatrix<U, Major>>;
        using Base::isReverseDiff;
    private:
        using ValueMatrix = device_obj<DenseMatrix<T, Major>>;
        using GradType = ScalarType::GradType;
        using GradMatrix = std::conditional<Order == 1, ValueMatrix, device_obj<DenseMatrix<GradType, Major>>>::type;

        ValueMatrix v;
        GradMatrix g;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(size_t row, size_t col, T init);
        device_obj(ValueMatrix v_, GradMatrix g_);
        template<Vector V>
        explicit(isReverseDiff()) device_obj(const V& vec) requires(!ReverseDiff<V>);
        template<Matrix M>
        explicit(isReverseDiff()) device_obj(const M& mat) requires(!ReverseDiff<M>);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        void zero_grad();

        using Base::resize;
        __host__ __device__ void resize(size_t row, size_t col);
        void reserve(size_t size);

        void zeros();
        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return v.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return v.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return v.getOrder(); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return v.empty(); }

        [[nodiscard]] __host__ __device__ auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto&& grads(this auto&&) noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t row, size_t col, auto& distribution);
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int Major>
    class Traits<device_obj<Diff<DenseMatrix<T, Major>, Mode, Order>>>
            : public Traits<DenseMatrix<T, Major>> {
    public:
        using ScalarType = device_obj<Diff<T, Mode, Order>>;
    };
}

#include "MatrixImpl/DiffDenseMatrixImpl.cuh"
