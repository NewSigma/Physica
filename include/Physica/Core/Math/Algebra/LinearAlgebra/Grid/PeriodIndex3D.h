/*
 * Copyright 2023 Weibo He.
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

#include <Physica/Core/Math/Transform/FFT.h>
#include "GridImpl/GridBase.h"

namespace Physica::Core {
    /**
     * Provide a algebra for Index1D and Index3D in periodic boundary condition
     */
    class PHYSICA_API PeriodIndex3D {
        using Index3D = typename GridBase::Index3D;
        using Index1D = typename Index3D::ElemType;

        Index3D index;
        Index3D dim;
    public:
        __host__ __device__ inline PeriodIndex3D(Index3D index_, Index3D dim_);
        __host__ __device__ inline PeriodIndex3D(Index1D index_, Index3D dim_);
        PeriodIndex3D(const PeriodIndex3D&) = default;
        PeriodIndex3D(PeriodIndex3D&&) noexcept = default;
        ~PeriodIndex3D() = default;
        /* Operators */
        __host__ __device__ inline PeriodIndex3D& operator=(PeriodIndex3D obj) noexcept;
        [[nodiscard]] __host__ __device__ operator Index3D() { return index; }
        [[nodiscard]] __host__ __device__ Index1D operator[](int i) { return index[i]; }
        [[nodiscard]] __host__ __device__ inline PeriodIndex3D operator+(const PeriodIndex3D& other) const;
        /* Operations */
        [[nodiscard]] __host__ __device__ inline Index1D toIndex1D() const;
        [[nodiscard]] __host__ __device__ inline bool isInReducedK() const;
        [[nodiscard]] __host__ __device__ inline Index3D toReducedK() const;
        __host__ __device__ inline void normalize();
        __host__ __device__ inline void swap(PeriodIndex3D& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Index3D getDim() const noexcept { return dim; }
    private:
        __host__ __device__ inline static Index3D toIndex3D(Index1D index, Index3D dim);
    };

    __host__ __device__ inline PeriodIndex3D::PeriodIndex3D(Index3D index_, Index3D dim_) : index(index_), dim(dim_) {
        for (int i = 0; i < 3; ++i)
            assert(dim[i] > 0 && "[Error]: Zero dim is not allowed");
    }

    __host__ __device__ inline PeriodIndex3D::PeriodIndex3D(Index1D index_, Index3D dim_) : PeriodIndex3D(toIndex3D(index_, dim_), dim_) {}

    __host__ __device__ inline PeriodIndex3D& PeriodIndex3D::operator=(PeriodIndex3D obj) noexcept {
        swap(obj);
        return *this;
    }

    __host__ __device__ inline PeriodIndex3D PeriodIndex3D::operator+(const PeriodIndex3D& other) const {
        assert(dim == other.dim && "[Error]: Inconsistent dimention");
        Index3D result{};
        for (int i = 0; i < 3; ++i)
            result[i] = (index[i] + other.index[i]) % dim[i];
        return PeriodIndex3D(result, dim);
    }

    __host__ __device__ inline typename PeriodIndex3D::Index1D PeriodIndex3D::toIndex1D() const {
        return (index[0] * dim[1] + index[1]) * dim[2] + index[2];
    }

    __host__ __device__ inline bool PeriodIndex3D::isInReducedK() const {
        const size_t kSpaceDimZ = FFT<Real<>>::rSizeToKSize(dim[2]);
        return index[2] < kSpaceDimZ;
    }

    __host__ __device__ inline typename PeriodIndex3D::Index3D PeriodIndex3D::toReducedK() const {
        Index3D result = index;
        if (!isInReducedK()) {
            result[0] = (dim[0] - result[0]) % dim[0];
            result[1] = (dim[1] - result[1]) % dim[1];
            result[2] = dim[2] - result[2];
        }
        return result;
    }

    __host__ __device__ inline void PeriodIndex3D::normalize() {
        for (int i = 0; i < 3; ++i)
            index[i] %= dim[i];
    }

    __host__ __device__ inline void PeriodIndex3D::swap(PeriodIndex3D& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        index.swap(obj.index);
        dim.swap(obj.dim);
    }

    __host__ __device__ inline typename PeriodIndex3D::Index3D PeriodIndex3D::toIndex3D(Index1D index, Index3D dim) {
        Index3D result{};
        result[0] = index / (dim[1] * dim[2]);
        const Index1D temp = index % (dim[1] * dim[2]);
        result[1] = temp / dim[2];
        result[2] = temp % dim[2];
        return result;
    }
}
