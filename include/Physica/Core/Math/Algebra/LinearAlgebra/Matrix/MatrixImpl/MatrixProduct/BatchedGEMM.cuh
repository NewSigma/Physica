/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstdint>
#include <utility>
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "BatchedGEMM.h"

namespace Physica {
    template<Scalar T>
    class device_obj<BatchedGEMM<T>> {
        using host_obj = BatchedGEMM<T>;
        using This = device_obj<host_obj>;
    public:
        using PointerTable = host_obj::PointerTable;
    private:
        PointerTable ptrs;
        host_obj::Shape shape;
        device_obj<PointerTable> d_ptrs;
        bool sync = false;
    public:
        device_obj() = default;
        device_obj(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, size_t batchSize);
        device_obj(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, PointerTable ptrs);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void batch(size_t index, Array<T*, 3> handles) noexcept;

        void compute(T alpha = T(1), T beta = T(0));
        void compute(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha = T(1), T beta = T(0));
        void compute_cublas(T alpha, T beta);
        void compute_cublas(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha, T beta);

        void reshape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c);
        void resize(size_t batchSize);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointers() const noexcept { return ptrs; }
        [[nodiscard]] const auto& getShape() const noexcept { return shape; }
        [[nodiscard]] size_t getBatchSize() const noexcept { return ptrs.getCol(); }
    };

    template<Scalar T>
    device_obj<BatchedGEMM<T>>::device_obj(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, size_t batchSize) {
        reshape(a, b, c);
        resize(batchSize);
    }

    template<Scalar T>
    device_obj<BatchedGEMM<T>>::device_obj(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, PointerTable ptrs_)
            : device_obj(a, b, c, ptrs_.getCol()) {
        for (size_t i = 0; i < ptrs_.getCol(); ++i)
            batch(i, {ptrs_[0, i], ptrs_[1, i], ptrs_[2, i]});
    }

    template<Scalar T>
    auto device_obj<BatchedGEMM<T>>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::batch(size_t index, Array<T*, 3> handles) noexcept {
        ptrs[0, index] = handles[shape.transC ? 1 : 0];
        ptrs[1, index] = handles[shape.transC ? 0 : 1];
        ptrs[2, index] = handles[2];
        sync = false;
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::compute(T alpha, T beta) {
        compute_cublas(alpha, beta);
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::compute(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha, T beta) {
        assert(getBatchSize() != 0 && "[Error]: Batch size must not be zero");
        compute_cublas(strideA, strideB, strideC, alpha, beta);
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::compute_cublas(T alpha_, T beta_) {
        using Tm = decltype(std::declval<T>().toCUDA());
        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        const Tm alpha = alpha_.toCUDA();
        const Tm beta = beta_.toCUDA();
        const auto opA = shape.transA ? CUBLAS_OP_T : CUBLAS_OP_N;
        const auto opB = shape.transB ? CUBLAS_OP_T : CUBLAS_OP_N;
        const int64_t m = int64_t(shape.row);
        const int64_t n = int64_t(shape.col);
        const int64_t k = int64_t(shape.inner);
        const int64_t lda = int64_t(shape.majorStrideA);
        const int64_t ldb = int64_t(shape.majorStrideB);
        const int64_t ldc = int64_t(shape.majorStrideC);
        const int64_t batchCount = int64_t(getBatchSize());
        if (!sync) {
            ptrs.toDeviceAsync(d_ptrs);
            sync = true;
        }
        auto* aArray = d_ptrs.data_ptr(0, 0);
        auto* bArray = d_ptrs.data_ptr(1, 0);
        auto* cArray = d_ptrs.data_ptr(2, 0);
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                check(cublasCgemmBatched_64(ctx, opA, opB, m, n, k, &alpha, (const cuComplex* const*)aArray, lda, (const cuComplex* const*)bArray, ldb, &beta, (cuComplex* const*)cArray, ldc, batchCount));
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                check(cublasZgemmBatched_64(ctx, opA, opB, m, n, k, &alpha, (const cuDoubleComplex* const*)aArray, lda, (const cuDoubleComplex* const*)bArray, ldb, &beta, (cuDoubleComplex* const*)cArray, ldc, batchCount));
            }
        }
        else {
            if constexpr (T::Prec == Float16)
                check(cublasHgemmBatched_64(ctx, opA, opB, m, n, k, &alpha, (const __half* const*)aArray, lda, (const __half* const*)bArray, ldb, &beta, (__half* const*)cArray, ldc, batchCount));
            else if constexpr (T::Prec == Float32)
                check(cublasSgemmBatched_64(ctx, opA, opB, m, n, k, &alpha, (const float* const*)aArray, lda, (const float* const*)bArray, ldb, &beta, (float* const*)cArray, ldc, batchCount));
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                check(cublasDgemmBatched_64(ctx, opA, opB, m, n, k, &alpha, (const double* const*)aArray, lda, (const double* const*)bArray, ldb, &beta, (double* const*)cArray, ldc, batchCount));
            }
        }
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::compute_cublas(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha_, T beta_) {
        using Tm = decltype(std::declval<T>().toCUDA());
        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        const Tm alpha = alpha_.toCUDA();
        const Tm beta = beta_.toCUDA();
        const auto opA = shape.transA ? CUBLAS_OP_T : CUBLAS_OP_N;
        const auto opB = shape.transB ? CUBLAS_OP_T : CUBLAS_OP_N;
        const int64_t m = int64_t(shape.row);
        const int64_t n = int64_t(shape.col);
        const int64_t k = int64_t(shape.inner);
        const int64_t lda = int64_t(shape.majorStrideA);
        const int64_t ldb = int64_t(shape.majorStrideB);
        const int64_t ldc = int64_t(shape.majorStrideC);
        const int64_t batchCount = int64_t(getBatchSize());
        const auto* aData = reinterpret_cast<const Tm*>(ptrs[0, 0]);
        const auto* bData = reinterpret_cast<const Tm*>(ptrs[1, 0]);
        auto* cData = reinterpret_cast<Tm*>(ptrs[2, 0]);
        const auto sa = static_cast<long long>(shape.transC ? strideB : strideA);
        const auto sb = static_cast<long long>(shape.transC ? strideA : strideB);
        const auto sc = static_cast<long long>(strideC);
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                check(cublasCgemmStridedBatched_64(ctx, opA, opB, m, n, k, (const cuComplex*)&alpha, (const cuComplex*)aData, lda, sa, (const cuComplex*)bData, ldb, sb, (const cuComplex*)&beta, (cuComplex*)cData, ldc, sc, batchCount));
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                check(cublasZgemmStridedBatched_64(ctx, opA, opB, m, n, k, (const cuDoubleComplex*)&alpha, (const cuDoubleComplex*)aData, lda, sa, (const cuDoubleComplex*)bData, ldb, sb, (const cuDoubleComplex*)&beta, (cuDoubleComplex*)cData, ldc, sc, batchCount));
            }
        }
        else {
            if constexpr (T::Prec == Float32)
                check(cublasSgemmStridedBatched_64(ctx, opA, opB, m, n, k, (const float*)&alpha, (const float*)aData, lda, sa, (const float*)bData, ldb, sb, (const float*)&beta, (float*)cData, ldc, sc, batchCount));
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                check(cublasDgemmStridedBatched_64(ctx, opA, opB, m, n, k, (const double*)&alpha, (const double*)aData, lda, sa, (const double*)bData, ldb, sb, (const double*)&beta, (double*)cData, ldc, sc, batchCount));
            }
        }
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::reshape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c) {
        shape = host_obj::makeShape(a, b, c);
    }
    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::resize(size_t batchSize) {
        assert(batchSize != 0 && "[Error]: Batch size must not be zero");
        if (getBatchSize() != batchSize) {
            ptrs.resize(3, batchSize);
            d_ptrs.resize(3, batchSize);
            sync = false;
        }
    }

    template<Scalar T>
    void device_obj<BatchedGEMM<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(shape, obj.shape);
        ptrs.swap(obj.ptrs);
        d_ptrs.swap(obj.d_ptrs);
        std::swap(sync, obj.sync);
    }

}
