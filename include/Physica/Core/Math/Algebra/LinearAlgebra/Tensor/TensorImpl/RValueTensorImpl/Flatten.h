/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"
#include "../RValueTensor.h"

namespace Physica {
    template<Tensor T>
    class Flatten<T> : public RValueVector<Flatten<T>> {
        using This = Flatten<T>;
        using Base = RValueVector<This>;

        decay_rvalue_t<T> tensor;
    public:
        Flatten(T&& tensor_) : tensor(std::forward<T>(tensor_)) {}
        Flatten(const This&) = default;
        Flatten(This&&) noexcept = default;
        ~Flatten() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(size_t index) const;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return tensor.getSize(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<T>::getSizeAtCompile(); }
    };

    template<Tensor T>
    decltype(auto) Flatten<T>::calc(size_t index) const {
        return tensor.calc(tensor.toIndexND(index));
    }
}
