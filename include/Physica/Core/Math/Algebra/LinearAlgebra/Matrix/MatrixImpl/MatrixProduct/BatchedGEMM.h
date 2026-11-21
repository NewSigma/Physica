/*
 * Copyright 2026 Weibo He.
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
        };
    private:
        PointerTable ptrs;
        Shape shape;
    public:
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
        void compute_base(T alpha, T beta) noexcept;
        void compute_mkl(T alpha, T beta) noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointers() const noexcept { return ptrs; }
        [[nodiscard]] size_t getBatchSize() const noexcept { return ptrs.getCol(); }
        [[nodiscard]] const auto& getShape() const noexcept { return shape; }
    private:
        [[nodiscard]] static auto makeShape(const Matrix auto& a, const Matrix auto& b, Matrix auto& c) -> Shape;
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    BatchedGEMM<T>::BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, size_t batchSize)
            : ptrs(3, batchSize)
            , shape(makeShape(a, b, c)) {
        assert(batchSize != 0 && "[Error]: Batch size must not be zero");
    }

    template<Scalar T>
    BatchedGEMM<T>::BatchedGEMM(const Matrix auto& a, const Matrix auto& b, Matrix auto& c, PointerTable ptrs_)
            : BatchedGEMM(a, b, c, ptrs_.getCol()) {
        ptrs = std::move(ptrs_);
    }

    template<Scalar T>
    auto BatchedGEMM<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void BatchedGEMM<T>::batch(size_t index, Array<T*, 3> handles) noexcept {
        ptrs[0, index] = handles[0];
        ptrs[1, index] = handles[1];
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
        constexpr bool TransA = instanceof<remove_device_obj_t<A>, Transpose>;
        constexpr bool TransB = instanceof<remove_device_obj_t<B>, Transpose>;
        static_assert(!instanceof<remove_device_obj_t<C>, Transpose>, "[Error]: The target must not be transposed");

        auto&& storageA = [&]() -> decltype(auto) {
            if constexpr (TransA)
                return std::forward<decltype(a)>(a).transpose();
            else
                return std::forward<decltype(a)>(a);
        }();
        auto&& storageB = [&]() -> decltype(auto) {
            if constexpr (TransB)
                return std::forward<decltype(b)>(b).transpose();
            else
                return std::forward<decltype(b)>(b);
        }();
        auto&& storageC = std::forward<decltype(c)>(c);

        Shape result{};
        result.row = a.getRow();
        result.inner = a.getCol();
        result.col = c.getCol();
        result.majorStrideA = storageA.getMajorStride();
        result.majorStrideB = storageB.getMajorStride();
        result.majorStrideC = storageC.getMajorStride();
        result.transA = TransA;
        result.transB = TransB;
        return result;
    }
}

#ifdef PHYSICA_MKL
    #include "BatchedGEMM_MKL.h"
#endif
