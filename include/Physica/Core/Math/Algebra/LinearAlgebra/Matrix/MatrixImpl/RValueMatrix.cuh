/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Parallel/ThreadBlock.cuh"
#include "RValueMatrixImpl/RMatrixBlock.cuh"
#include "RValueMatrix.h"

namespace Physica {
    template<class Derived>
    class device_obj<RValueMatrix<Derived>> : public CRTPBase<device_obj<RValueMatrix<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is unnecessary");
        using host_obj = RValueMatrix<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<device_obj<Derived>>;
    public:
        using ScalarType = TraitsType::ScalarType;
        constexpr static int Major = TraitsType::Major;
        constexpr static size_t RowAtCompile = TraitsType::RowAtCompile;
        constexpr static size_t ColAtCompile = TraitsType::ColAtCompile;
        constexpr static size_t SizeAtCompile = TraitsType::SizeAtCompile;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isDiffable = ScalarType::isDiffable;
        constexpr static bool isComplex = ScalarType::isComplex;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] __host__ __device__ auto operator*(this auto&&, Vector auto&& v) noexcept requires(RowAtCompile != 1);
        [[nodiscard]] __host__ __device__ auto operator*(this auto&&, Matrix auto&& m) noexcept;
        /* Operations */
        __host__ __device__ void assign(Matrix auto&& target) const;
        __device__ void assign(Matrix auto&& target, const ThreadBlock& block) const;
        __host__ __device__ void assign_add(Matrix auto&& target) const;
        __device__ void assign_add(Matrix auto&& target, const ThreadBlock& block) const;
        __host__ __device__ void assert_assign(const Matrix auto& source) const noexcept;

        [[nodiscard]] __device__ auto calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] __device__ auto calc_value(size_t row, size_t col) const { return Base::getDerived().calc_value(row, col); }
        [[nodiscard]] __device__ T calcFromMajorMinor(size_t major, size_t minor) const;
        void reverse(const Matrix auto& y, const Matrix auto& grad) const noexcept;

        __host__ __device__ void resize(const Matrix auto& m, auto&&... args);
        __host__ __device__ auto resize(size_t r, size_t c, auto&&... args);

        [[nodiscard]] __host__ __device__ auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] __host__ __device__ auto col(this auto&&, size_t c) noexcept;
        [[nodiscard]] __host__ __device__ auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] __host__ __device__ auto topRows(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRows(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ auto leftCols(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto rightCols(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ auto diag(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto triu(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto triu_unit(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto tril(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto tril_unit(this auto&&) noexcept;

        [[nodiscard]] __host__ __device__ auto sum_rows() const;
        [[nodiscard]] __host__ __device__ auto sum_cols() const;

        [[nodiscard]] __host__ __device__ auto inv(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto transpose(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto hermite() const noexcept;
        [[nodiscard]] __host__ __device__ auto flatten() const noexcept;

        [[nodiscard]] __host__ __device__ decltype(auto) reals(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto imags(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto squaredNorms(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto norms(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grads(this auto&&) noexcept;

        [[nodiscard]] __host__ __device__ KernelConfig makeKernelConfig() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getCol(); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getRow() * getCol(); }
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept { return MatrixMajor::getMaxMajor<device_obj<Derived>>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept { return MatrixMajor::getMaxMinor<device_obj<Derived>>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return Base::getDerived().getSize() == 0; }

        [[nodiscard]] __host__ __device__ bool isSquare() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept { return host_obj::isCompact(); }
        [[nodiscard]] __host__ __device__ static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept;
        [[nodiscard]] __host__ __device__ static size_t colFromMajorMinor(size_t major, size_t minor) noexcept;
        [[nodiscard]] __host__ __device__ static KernelConfig makeKernelConfig(size_t maxMajor, size_t maxMinor) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Scalar auto& source) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Matrix auto& source) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RValueMatrix<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.cuh"
#include "RValueMatrixImpl/Sum.cuh"
#include "RValueMatrixImpl/DiagVector.cuh"
#include "RValueMatrixImpl/Transpose.cuh"
#include "RValueMatrixImpl/Hermite.cuh"
#include "RValueMatrixImpl/Convert/MatrixConvert.cuh"
#include "MatrixProduct/GEMM.cuh"
#include "MatrixProduct/GEMV.cuh"
#include "MatrixProduct/GEVM.cuh"
#include "MatrixProduct/Kronecker.cuh"
#include "MatrixExpr.cuh"
