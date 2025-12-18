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
    template<Matrix M1, Matrix M2>
    class device_obj<GEMM<M1, M2>> : public device_obj<RValueMatrix<GEMM<M1, M2>>> {
        using host_obj = GEMM<M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tv;
    private:
        using Tm = ScalarType::MachineType;

        Physica::PlainStruct<const M1> mat1;
        Physica::PlainStruct<const M2> mat2;
    public:
        __host__ __device__ device_obj(const M1& mat1_, const M2& mat2_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void assign(Matrix auto& target) const;
        __host__ __device__ void assign_add(Matrix auto& target) const;

        [[nodiscard]] auto compute() const;
        [[nodiscard]] __device__ ScalarType calc(size_t, size_t) const { noImpl("GEMM.calc() is low performance and should be avoided"); }
        [[nodiscard]] __device__ Tv calc_value(size_t, size_t) const { noImpl("GEMM.calc_value() is low performance and should be avoided"); }

        void reverse(const Matrix auto& grad) const noexcept;

        __host__ __device__ auto values() const noexcept { return getLHS().values() * getRHS().values(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getRHS().getCol(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return mat1.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return mat2.getDerived(); }
    private:
        template<bool AssignAdd>
        void assign_impl_cublas(Matrix auto& target) const;
    };

    template<Matrix M1, Matrix M2>
    __host__ __device__ device_obj<GEMM<M1, M2>>::device_obj(const M1& mat1_, const M2& mat2_) : mat1(asStruct(mat1_)), mat2(asStruct(mat2_)) {
        assert(mat1_.getCol() == mat2_.getRow());
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ void device_obj<GEMM<M1, M2>>::assign(Matrix auto& target) const {
        if constexpr (IsHost())
            assign_impl_cublas<false>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ void device_obj<GEMM<M1, M2>>::assign_add(Matrix auto& target) const {
        if constexpr (IsHost())
            assign_impl_cublas<true>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix M1, Matrix M2>
    auto device_obj<GEMM<M1, M2>>::compute() const {
        using RetTy = decltype(std::declval<host_obj>().compute())::device_obj_type;
        return RetTy(*this);
    }

    template<Matrix M1, Matrix M2>
    void device_obj<GEMM<M1, M2>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        if constexpr (ReverseDiff<M1>)
            getLHS().reverse(grad * getRHS().values().transpose());
        if constexpr (ReverseDiff<M2>)
            getRHS().reverse(getLHS().values().transpose() * grad);
    }

    template<Matrix M1, Matrix M2>
    template<bool AssignAdd>
    void device_obj<GEMM<M1, M2>>::assign_impl_cublas(Matrix auto& target) const {
        constexpr bool IsDeviceMatrix = Traits<M1>::SizeAtCompile == Dynamic && Traits<M2>::SizeAtCompile == Dynamic;
        constexpr bool isTranspose1 = is_transpose<M1>::value;
        constexpr bool isTranspose2 = is_transpose<M2>::value;
        static_assert(IsDeviceMatrix, "[Error]: Fixed matrix is on host, pass it to device before calling cuBLAS");
        static_assert(MatrixOption::isColMatrix<M1>() && MatrixOption::isColMatrix<M2>(), "[Error]: cuBLAS uses column major");
        static_assert(!Diffable<This>);
        target.assert_assign(*this);

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

        const Tm* A{};
        size_t k{};
        if constexpr (isTranspose1) {
            k = getLHS().getRow();
            A = (Tm*)getLHS().getExpr().data();
        }
        else {
            k = getLHS().getCol();
            A = (Tm*)getLHS().data();
        }

        const Tm* B{};
        if constexpr (isTranspose2)
            B = (Tm*)getRHS().getExpr().data();
        else
            B = (Tm*)getRHS().data();

        Tm* C{};
        if constexpr (Diffable<decltype(target)>)
            C = (Tm*)target.data().value_ptr();
        else
            C = (Tm*)target.data();

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (ScalarType::Prec == Float16)
            check(cublasHgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else if constexpr (ScalarType::Prec == Float32)
            check(cublasSgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else {
            static_assert(ScalarType::Prec == Float64, "[Error]: Unknown ScalarType");
            check(cublasDgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        }
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] __host__ __device__ auto operator*(const M1& mat1, const M2& mat2) noexcept
            requires(((M1::ColAtCompile != 1 && M2::ColAtCompile != 1) || (M1::ColAtCompile == 1 && M2::ColAtCompile == 1)) && CUDA<M1> && CUDA<M2>) {
        return device_obj<GEMM<M1, M2>>(mat1, mat2);
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<device_obj<GEMM<M1, M2>>> : public Traits<GEMM<M1, M2>> {
    public:
        constexpr static bool FastAssign = true;
    };
}
