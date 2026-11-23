/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>
#include "Physica/Core/Utils/Container/Array2D.h"
#include "Physica/Core/Utils/NoImpl.h"

namespace Physica {
    template<Scalar T>
    class BatchedGEMM {
        using This = BatchedGEMM;
    public:
        using PointerTable = Array2D<T*, MatrixMajor::Row, 3>;
        struct Shape {
            size_t row = 0;
            size_t col = 0;
            size_t inner = 0;
            size_t majorStrideA = 0;
            size_t majorStrideB = 0;
            size_t majorStrideC = 0;
            bool transA = false;
            bool transB = false;
            bool transC = false;
        };
    private:
        PointerTable ptrs;
        Shape shape;
    public:
        BatchedGEMM() = default;
        BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, size_t batchSize);
        BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, PointerTable ptrs_);
        BatchedGEMM(const This&) = default;
        BatchedGEMM(This&&) noexcept = default;
        ~BatchedGEMM() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void batch(size_t index, Array<T*, 3> handles) noexcept;

        void compute(T alpha = T(1), T beta = T(0));
        void compute(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha = T(1), T beta = T(0));
        void compute_base(T alpha, T beta) noexcept;
        void compute_mkl(T alpha, T beta) noexcept;
        void compute_mkl(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha, T beta) noexcept;

        void reshape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c);
        void resize(size_t batchSize);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointers() const noexcept { return ptrs; }
        [[nodiscard]] size_t getBatchSize() const noexcept { return ptrs.getCol(); }
        [[nodiscard]] const auto& getShape() const noexcept { return shape; }
    private:
        [[nodiscard]] static auto makeShape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c) -> Shape;
        [[nodiscard]] bool isPackedStrided(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC) const noexcept;
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    BatchedGEMM<T>::BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, size_t batchSize) {
        reshape(a, b, c);
        resize(batchSize);
    }

    template<Scalar T>
    BatchedGEMM<T>::BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, PointerTable ptrs_)
            : BatchedGEMM(a, b, c, ptrs_.getCol()) {
        for (size_t i = 0; i < ptrs_.getCol(); ++i)
            batch(i, {ptrs_[0, i], ptrs_[1, i], ptrs_[2, i]});
    }

    template<Scalar T>
    auto BatchedGEMM<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void BatchedGEMM<T>::batch(size_t index, Array<T*, 3> handles) noexcept {
        ptrs[0, index] = handles[shape.transC ? 1 : 0];
        ptrs[1, index] = handles[shape.transC ? 0 : 1];
        ptrs[2, index] = handles[2];
    }

    template<Scalar T>
    void BatchedGEMM<T>::compute(T alpha, T beta) {
        if constexpr (HasMKL() && !T::isDiffable() && (T::Prec == Float32 || T::Prec == Float64))
            compute_mkl(alpha, beta);
        else
            compute_base(alpha, beta);
    }

    template<Scalar T>
    void BatchedGEMM<T>::compute_base(T, T) noexcept {
        noImpl("No host batched GEMM support");
    }

    template<Scalar T>
    void BatchedGEMM<T>::compute(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC, T alpha, T beta) {
        const size_t batchSize = getBatchSize();
        assert(batchSize != 0 && "[Error]: Batch size must not be zero");
        if constexpr (HasMKL() && !T::isDiffable() && (T::Prec == Float32 || T::Prec == Float64)) {
            if (!isPackedStrided(strideA, strideB, strideC)) {
                // MKL does not support; fallback to pointer table
                T* const first = ptrs[0, 0];
                T* const second = ptrs[1, 0];
                T* const baseA = shape.transC ? second : first;
                T* const baseB = shape.transC ? first : second;
                T* const baseC = ptrs[2, 0];
                for (size_t index = 0; index < batchSize; ++index) {
                    const auto offset = std::ptrdiff_t(index);
                    batch(index, {baseA + offset * strideA, baseB + offset * strideB, baseC + offset * strideC});
                }
                compute_mkl(alpha, beta);
                return;
            }
            compute_mkl(strideA, strideB, strideC, alpha, beta);
        }
        else
            noImpl("No host strided batched GEMM support");
    }

    template<Scalar T>
    void BatchedGEMM<T>::reshape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c) {
        shape = makeShape(a, b, c);
    }

    template<Scalar T>
    void BatchedGEMM<T>::resize(size_t batchSize) {
        assert(batchSize != 0 && "[Error]: Batch size must not be zero");
        if (getBatchSize() != batchSize)
            ptrs.resize(3, batchSize);
    }

    template<Scalar T>
    void BatchedGEMM<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(shape, obj.shape);
        ptrs.swap(obj.ptrs);
    }

    template<Scalar T>
    auto BatchedGEMM<T>::makeShape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c) -> Shape {
        using A = std::remove_cvref_t<decltype(a)>;
        using B = std::remove_cvref_t<decltype(b)>;
        using C = std::remove_cvref_t<decltype(c)>;
        static_assert(A::isStrided(), "[Error]: The operands must be strided matrices");
        static_assert(B::isStrided(), "[Error]: The operands must be strided matrices");
        static_assert(C::isStrided(), "[Error]: The operands must be strided matrices");
        static_assert(A::getMinorStrideAtCompile() == 1
                   && B::getMinorStrideAtCompile() == 1
                   && C::getMinorStrideAtCompile() == 1,
                      "[Error]: BLAS expects a contiguous minor direction");
        constexpr bool IsTransposeA = instanceof<remove_device_obj_t<A>, Transpose>;
        constexpr bool IsTransposeB = instanceof<remove_device_obj_t<B>, Transpose>;
        constexpr bool TransC = MatrixMajor::isRowMatrix<remove_device_obj_t<C>>();

        auto&& storageA = [&]() -> decltype(auto) {
            if constexpr (IsTransposeA)
                return std::forward<decltype(a)>(a).transpose();
            else
                return std::forward<decltype(a)>(a);
        }();
        auto&& storageB = [&]() -> decltype(auto) {
            if constexpr (IsTransposeB)
                return std::forward<decltype(b)>(b).transpose();
            else
                return std::forward<decltype(b)>(b);
        }();

        Shape result{};
        if constexpr (TransC) {
            result.row = c.getCol();
            result.col = c.getRow();
            result.inner = a.getCol();
            result.transA = MatrixMajor::isColMatrix<B>();
            result.transB = MatrixMajor::isColMatrix<A>();
            result.majorStrideA = storageB.getMajorStride();
            result.majorStrideB = storageA.getMajorStride();
        }
        else {
            result.row = a.getRow();
            result.inner = a.getCol();
            result.col = c.getCol();
            result.transA = MatrixMajor::isRowMatrix<A>();
            result.transB = MatrixMajor::isRowMatrix<B>();
            result.majorStrideA = storageA.getMajorStride();
            result.majorStrideB = storageB.getMajorStride();
        }
        result.majorStrideC = c.getMajorStride();
        result.transC = TransC;
        return result;
    }

    template<Scalar T>
    auto BatchedGEMM<T>::isPackedStrided(std::ptrdiff_t strideA, std::ptrdiff_t strideB, std::ptrdiff_t strideC) const noexcept -> bool {
        const std::ptrdiff_t stridea = shape.transC ? strideB : strideA;
        const std::ptrdiff_t strideb = shape.transC ? strideA : strideB;
        const std::ptrdiff_t colsA = shape.transA ? shape.row : shape.inner;
        const std::ptrdiff_t colsB = shape.transB ? shape.inner : shape.col;
        return stridea >= std::ptrdiff_t(shape.majorStrideA * colsA)
            && strideb >= std::ptrdiff_t(shape.majorStrideB * colsB)
            && strideC >= std::ptrdiff_t(shape.majorStrideC * shape.col);
    }
}

#ifdef PHYSICA_MKL
    #include "BatchedGEMM_MKL.h"
#endif
