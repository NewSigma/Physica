/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/ContinuousMatrix.cuh"
#include "DenseMatrix.cuh"
#include "DiffDenseMatrix.h"

namespace Physica {
    template<Scalar T, int Order, int Option>
    class device_obj<DenseMatrix<Diff<T, DiffMode::Reverse, Order>, Option>> : public device_obj<ContinuousMatrix<DenseMatrix<Diff<T, DiffMode::Reverse, Order>, Option>>> {
        static_assert(!T::isDiffable, "[Error]: Nested Diff<> is not allowed");
        using host_obj = DenseMatrix<Diff<T, DiffMode::Reverse, Order>, Option>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousMatrix<host_obj>>;
    public:
        using ScalarType = Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using ValueMatrix = device_obj<DenseMatrix<T, Option>>;
        using GradType = ScalarType::GradType;
        using GradMatrix = std::conditional<Order == 1, ValueMatrix, device_obj<DenseMatrix<GradType, Option>>>::type;

        ValueMatrix v;
        GradMatrix g;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(size_t row, size_t col, T init);
        device_obj(ValueMatrix v_, GradMatrix g_);
        template<Vector V>
        explicit(isReverseDiff) device_obj(const V& vec) requires(!ReverseDiff<V>);
        template<Matrix M>
        explicit(isReverseDiff) device_obj(const M& mat) requires(!ReverseDiff<M>);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator();
        /* Operations */
        template<RNG R>
        inline void random_uniform();
        template<RNG R>
        inline void random_normal();
        template<RNG R, class Distribution>
        inline void random_any(Distribution& dist);

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return v.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return v.getRow(); }

        [[nodiscard]] __host__ __device__ const auto& values() const noexcept { return v; }
        [[nodiscard]] __host__ __device__ auto& values() noexcept { return v; }
        [[nodiscard]] __host__ __device__ const auto& grads() const noexcept { return g; }
        [[nodiscard]] __host__ __device__ auto& grads() noexcept { return g; }
        /* Static members */
        template<RNG R>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] inline static This random_normal(size_t row, size_t col);
        template<RNG R, class Distribution>
        [[nodiscard]] inline static This random_any(size_t row, size_t col, Distribution& dist);
    private:
    };
}

namespace Physica {
    template<Scalar T, int Order, int Option>
    class Traits<device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>>
            : public Traits<DenseMatrix<T, Option>> {
    public:
        using ScalarType = device_obj<Diff<T, DiffMode::Reverse, Order>>;
    };
}

#include "MatrixImpl/DiffDenseMatrixImpl.cuh"
