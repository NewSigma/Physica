/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../LValueTensor.h"

namespace Physica {
    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    class TensorFiber<X, Dim> : public LValueVector<TensorFiber<X, Dim>> {
        using This = TensorFiber<X, Dim>;
        using Base = LValueVector<TensorFiber<X, Dim>>;
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
        [[nodiscard]] auto data_ptr(this auto&& self, size_t i) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getShapeAtCompile()[Dim]; }
    };

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    TensorFiber<X, Dim>::TensorFiber(X&& tensor, IndexVar auto... indices) : tensor(std::forward<X>(tensor)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    void TensorFiber<X, Dim>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Tensor X, int Dim> requires(std::remove_cvref_t<X>::isLValueTensor() && !std::remove_cvref_t<X>::isStrided())
    auto TensorFiber<X, Dim>::data_ptr(this auto&& self, size_t i) noexcept {
        auto idx = self.index;
        idx[Dim] = i;
        return self.tensor.data_ptr(idx);
    }
}
