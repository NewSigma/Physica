/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.h"
#include "../CompactTensor.h"

namespace Physica {
    template<Tensor T> requires(std::remove_cvref_t<T>::isCompact())
    class Flatten<T> : public CompactVector<Flatten<T>> {
        using This = Flatten<T>;

        decay_rvalue_t<T> tensor;
    public:
        using Base = CompactVector<This>;
    public:
        Flatten(T&& tensor_) : tensor(std::forward<T>(tensor_)) {}
        Flatten(const This&) = default;
        Flatten(This&&) noexcept = default;
        ~Flatten() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return tensor.getSize(); }
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<T>::getSizeAtCompile(); }
    };

    template<Tensor T> requires(std::remove_cvref_t<T>::isCompact())
    auto Flatten<T>::data_handle(this auto&& self) noexcept {
        return self.tensor.data_handle();
    }
}
