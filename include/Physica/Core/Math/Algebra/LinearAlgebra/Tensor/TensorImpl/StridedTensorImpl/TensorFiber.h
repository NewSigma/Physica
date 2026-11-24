/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../CompactTensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.h"

namespace Physica {
    namespace Internal {
        template<class X, int Dim>
        using TensorFiberBase = std::conditional_t<std::remove_cvref_t<X>::getStrideAtCompile()[Dim] == 1,
                                                   CompactVector<TensorFiber<X, Dim>>,
                                                   StridedVector<TensorFiber<X, Dim>>>;
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    class TensorFiber<X, Dim> : public Internal::TensorFiberBase<X, Dim> {
        using This = TensorFiber<X, Dim>;
        using Base = Internal::TensorFiberBase<X, Dim>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(Dim < std::remove_cvref_t<X>::ndim());
    private:
        decay_rvalue_t<X> tensor;
        IndexType index;
    public:
        TensorFiber(X&& tensor, IndexVar auto... indices);
        TensorFiber(const This&) = default;
        TensorFiber(This&&) noexcept = default;
        ~TensorFiber() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t length);
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return tensor.dim(Dim); }
        [[nodiscard]] constexpr size_t getStride() const noexcept;
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getShapeAtCompile()[Dim]; }
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept;
    };

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    TensorFiber<X, Dim>::TensorFiber(X&& tensor, IndexVar auto... indices) : tensor(std::forward<X>(tensor)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(static_cast<size_t>(indices) < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    void TensorFiber<X, Dim>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    constexpr size_t TensorFiber<X, Dim>::getStride() const noexcept {
        return tensor.getStride(Dim);
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    auto TensorFiber<X, Dim>::data_handle(this auto&& self) noexcept {
        auto idx = self.index;
        idx[Dim] = 0;
        return self.tensor.data_ptr(idx);
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isStrided())
    __host__ __device__ consteval size_t TensorFiber<X, Dim>::getStrideAtCompile() noexcept {
        return std::remove_cvref_t<X>::getStrideAtCompile()[Dim];
    }
}
