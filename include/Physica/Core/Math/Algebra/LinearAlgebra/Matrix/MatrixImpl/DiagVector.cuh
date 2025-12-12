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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.cuh"
#include "RValueMatrix.cuh"

namespace Physica {
    template<Matrix M, bool isLValueMatrix>
    class device_obj<DiagVector<M, isLValueMatrix>> : public device_obj<RValueVector<DiagVector<M, isLValueMatrix>>> {
        using host_obj = DiagVector<M, isLValueMatrix>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<const device_obj<M>> mat;
    public:
        explicit device_obj(const device_obj<M>& mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __host__ __device__ const auto& getExpr() const noexcept { return mat.getDerived(); }
        [[nodiscard]] __device__ T calc(size_t index) const { return getExpr().calc(index, index); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getRow(); }
    };
}

namespace Physica {
    template<Matrix M, bool isLValueMatrix>
    class Traits<device_obj<DiagVector<M, isLValueMatrix>>> : public Traits<DiagVector<M, isLValueMatrix>> {};
}
