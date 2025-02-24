/*
 * Copyright 2023-2024 Weibo He.
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

        T& tensor;
    public:
        using Base = LValueVector<FlattenL<T>>;
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        FlattenL(T& tensor_) : tensor(tensor_) {}
        FlattenL(const This&) = delete;
        FlattenL(This&&) noexcept = delete;
        ~FlattenL() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return tensor.getSize(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    };

    template<Tensor T>
    __host__ __device__ auto FlattenL<T>::data_ptr(size_t index) -> PtrTy {
        return tensor.data_ptr(tensor.toIndexND(index));
    }

    template<Tensor T>
    __host__ __device__ auto FlattenL<T>::data_ptr(size_t index) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<Tensor T>
    class Traits<FlattenL<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
