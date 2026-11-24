/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"
#include "../RValueTensor.h"

namespace Physica {
    template<Tensor T> requires(std::remove_cvref_t<T>::isLValueTensor() && !std::remove_cvref_t<T>::isCompact())
    class Flatten<T> : public LValueVector<Flatten<T>> {
        using This = Flatten<T>;

        decay_rvalue_t<T> tensor;
    public:
        using Base = LValueVector<Flatten<T>>;
        using typename Base::ScalarType;
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
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<T>::getSizeAtCompile(); }
    };

    template<Tensor T> requires(std::remove_cvref_t<T>::isLValueTensor() && !std::remove_cvref_t<T>::isCompact())
    auto Flatten<T>::data_ptr(this auto&& self, size_t index) noexcept {
        return self.tensor.data_ptr(self.tensor.toIndexND(index));
    }
}

namespace Physica {
    template<Tensor T>
    class Traits<Flatten<T>> {
    public:
        using ScalarType = std::remove_cvref_t<T>::ScalarType;
    };
}
