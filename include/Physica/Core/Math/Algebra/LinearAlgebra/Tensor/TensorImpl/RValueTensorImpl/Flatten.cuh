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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.cuh"
#include "../RValueTensor.cuh"
#include "Flatten.h"

namespace Physica {
    template<Tensor X>
    class device_obj<Flatten<X>> : public device_obj<RValueVector<Flatten<X>>> {
        using host_obj = Flatten<X>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj_t<X>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<X>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return tensor.getDerived().getSize(); }
    };

    template<Tensor X>
    __device__ auto device_obj<Flatten<X>>::calc(size_t index, [[maybe_unused]] instanceof_x<ThreadBlock> auto block) const -> T {
        assert(index < getLength());
        return tensor.getDerived().calc(tensor.getDerived().toIndexND(index));
    }
}

namespace Physica {
    template<Tensor X>
    class Traits<device_obj<Flatten<X>>> : public Traits<Flatten<X>> {};
}
