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

#include "Physica/PlainStruct.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "GEMM.h"

namespace Physica {
    template<Matrix T1, Matrix T2>
    class device_obj<MatrixProduct<T1, T2>> : public device_obj<RValueMatrix<MatrixProduct<T1, T2>>> {
        using host_obj = MatrixProduct<T1, T2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Tv;
        using Base::isReverseDiff;
        using DefaultType = host_obj::DefaultType::device_obj_type;
    private:
        using Tm = ScalarType::MachineType;

        Physica::PlainStruct<const T1> mat1;
        Physica::PlainStruct<const T2> mat2;
    public:
        __host__ __device__ device_obj(const T1& mat1_, const T2& mat2_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        __host__ __device__ void assign(M& target) const requires(CUDA<M>);
        template<Matrix M>
        __host__ __device__ void assign_add(M& target) const requires(CUDA<M>);

        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }
        [[nodiscard]] __device__ ScalarType calc(size_t, size_t) const { noImpl("GEMM.calc() is low performance and should be avoided"); }
        [[nodiscard]] __device__ Tv calc_value(size_t, size_t) const { noImpl("GEMM.calc_value() is low performance and should be avoided"); }

        template<Matrix M>
        void reverse(const M& grad) const noexcept requires(isReverseDiff);

        __host__ __device__ auto values() const noexcept { return getLHS().values() * getRHS().values(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getRHS().getCol(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return mat1.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return mat2.getDerived(); }
    private:
        template<Matrix M, bool AssignAdd>
        void assign_impl_cublas(M& target) const requires(CUDA<M>);
    };

    template<Matrix T1, Matrix T2>
    __host__ __device__ device_obj<MatrixProduct<T1, T2>>::device_obj(const T1& mat1_, const T2& mat2_) : mat1(asStruct(mat1_)), mat2(asStruct(mat2_)) {
        assert(mat1_.getCol() == mat2_.getRow());
    }

    template<Matrix T1, Matrix T2>
    template<Matrix M>
    __host__ __device__ void device_obj<MatrixProduct<T1, T2>>::assign(M& target) const requires(CUDA<M>) {
        if constexpr (IsHost())
            assign_impl_cublas<M, false>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix T1, Matrix T2>
    template<Matrix M>
    __host__ __device__ void device_obj<MatrixProduct<T1, T2>>::assign_add(M& target) const requires(CUDA<M>) {
        if constexpr (IsHost())
            assign_impl_cublas<M, true>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix T1, Matrix T2>
    template<Matrix M>
    void device_obj<MatrixProduct<T1, T2>>::reverse(const M& grad) const noexcept requires(isReverseDiff) {
        if constexpr (ReverseDiff<T1>)
            getLHS().reverse(grad * getRHS().values().transpose());
        if constexpr (ReverseDiff<T2>)
            getRHS().reverse(getLHS().values().transpose() * grad);
    }

    template<Matrix T1, Matrix T2>
    template<Matrix M, bool AssignAdd>
    void device_obj<MatrixProduct<T1, T2>>::assign_impl_cublas(M& target) const requires(CUDA<M>) {
        using U1 = remove_transpose<T1>::Type;
        using U2 = remove_transpose<T2>::Type;
        constexpr bool IsDeviceMatrix = Traits<T1>::SizeAtCompile == Dynamic && Traits<T2>::SizeAtCompile == Dynamic;
        constexpr bool isTranspose1 = is_transpose<T1>::value;
        constexpr bool isTranspose2 = is_transpose<T2>::value;
        static_assert(IsDeviceMatrix, "[Error]: Fixed matrix is on host, pass it to device before calling cuBLAS");
        static_assert(MatrixOption::isColMatrix<U1>() && MatrixOption::isColMatrix<U2>(), "[Error]: cuBLAS uses column major");
        static_assert(MatrixOption::isElementMatrix<U1>() && MatrixOption::isElementMatrix<U2>(), "[Error]: cuBLAS need element storage");
        static_assert(!Diffable<This>);

        using Tm = ScalarType::MachineType;
        const auto op1 = isTranspose1 ? CUBLAS_OP_T : CUBLAS_OP_N;
        const auto op2 = isTranspose2 ? CUBLAS_OP_T : CUBLAS_OP_N;
        const size_t r = getRow();
        const size_t c = getCol();
        const Tm alpha = 1;
        const Tm* pAlpha = &alpha;

        Tm beta;
        if constexpr (AssignAdd)
            beta = 1;
        else
            beta = 0;
        const Tm* pBeta = &beta;

        const Tm* A;
        size_t k;
        if constexpr (isTranspose1) {
            k = getLHS().getRow();
            A = (Tm*)getLHS().getExpr().data();
        }
        else {
            k = getLHS().getCol();
            A = (Tm*)getLHS().data();
        }

        const Tm* B;
        if constexpr (isTranspose2)
            B = (Tm*)getRHS().getExpr().data();
        else
            B = (Tm*)getRHS().data();

        Tm* C;
        if constexpr (Diffable<M>)
            C = (Tm*)target.data().value_ptr();
        else
            C = (Tm*)target.data();

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (ScalarType::Option == Float16)
            check(cublasHgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else if constexpr (ScalarType::Option == Float32)
            check(cublasSgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else {
            static_assert(ScalarType::Option == Float64, "[Error]: Unknown ScalarType");
            check(cublasDgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        }
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] __host__ __device__ inline auto operator*(const T1& mat1, const T2& mat2) noexcept
            requires(((T1::ColAtCompile != 1 && T2::ColAtCompile != 1) || (T1::ColAtCompile == 1 && T2::ColAtCompile == 1)) && CUDA<T1> && CUDA<T2>) {
        return device_obj<MatrixProduct<T1, T2>>(mat1, mat2);
    }
}

namespace Physica {
    template<Matrix T1, Matrix T2>
    class Traits<device_obj<MatrixProduct<T1, T2>>> : public Traits<MatrixProduct<T1, T2>> {
    public:
        constexpr static bool FastAssign = true;
    };
}
