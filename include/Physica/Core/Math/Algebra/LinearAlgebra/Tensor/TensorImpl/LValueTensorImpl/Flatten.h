/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"
#include "../RValueTensor.h"

namespace Physica {
    template<Tensor T>
    class FlattenL<T> : public LValueVector<FlattenL<T>> {
        using This = FlattenL<T>;

        decay_rvalue_t<T> tensor;
    public:
        using Base = LValueVector<FlattenL<T>>;
        using typename Base::ScalarType;
    public:
        FlattenL(T&& tensor_) : tensor(std::forward<T>(tensor_)) {}
        FlattenL(const This&) = default;
        FlattenL(This&&) noexcept = default;
        ~FlattenL() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return tensor.getSize(); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Tensor T>
    auto FlattenL<T>::data_ptr(this auto&& self, size_t index) noexcept {
        return self.tensor.data_ptr(self.tensor.toIndexND(index));
    }
}

namespace Physica {
    template<Tensor T>
    class Traits<FlattenL<T>> {
    public:
        using ScalarType = std::remove_cvref_t<T>::ScalarType;
    };
}
