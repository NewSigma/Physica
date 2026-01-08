/*
 * Copyright 2025 Weibo He.
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

#include "../LValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<DiagVectorL<M>> : public device_obj<LValueVector<DiagVectorL<M>>> {
        using host_obj = DiagVectorL<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
    public:
        explicit device_obj(const device_obj<M>& mat) : mat(asStruct(mat)) {}
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] __host__ __device__ const auto& getExpr() const noexcept { return mat.getDerived(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getRow(); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept { return self.getExpr().data_ptr(index, index); }
    };
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<DiagVectorL<M>>> : public Traits<DiagVectorL<M>> {};
}
