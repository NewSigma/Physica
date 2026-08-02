/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica {
    /**
     * Provide a algebra for Index1D and Index3D in periodic boundary condition
     */
    class PHYSICA_API PeriodIndex3D {
        using This = PeriodIndex3D;
        using Index1D = Index3D::ElemType;

        Index3D index;
        Index3D shape;
    public:
        __host__ __device__ PeriodIndex3D(Index3D index_, Index3D shape_);
        __host__ __device__ PeriodIndex3D(Index1D index_, Index3D shape_);
        PeriodIndex3D(const This&) = default;
        PeriodIndex3D(This&&) noexcept = default;
        ~PeriodIndex3D() = default;
        /* Operators */
        __host__ __device__ This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ operator Index3D() { return index; }
        [[nodiscard]] __host__ __device__ Index1D operator[](int i) { return index[i]; }
        [[nodiscard]] __host__ __device__ This operator+(const This& other) const;
        /* Operations */
        [[nodiscard]] __host__ __device__ Index1D toIndex1D() const;
        [[nodiscard]] __host__ __device__ bool isInReducedK() const;
        [[nodiscard]] __host__ __device__ Index3D toReducedK() const;
        __host__ __device__ void normalize();
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Index3D getShape() const noexcept { return shape; }
    private:
        __host__ __device__ static Index3D toIndex3D(Index1D index, Index3D shape);
    };

    __host__ __device__ inline PeriodIndex3D::PeriodIndex3D(Index3D index_, Index3D shape_) : index(index_), shape(shape_) {
        for (int i = 0; i < 3; ++i)
            assert(shape[i] > 0 && "[Error]: Zero shape is not allowed");
    }

    __host__ __device__ inline PeriodIndex3D::PeriodIndex3D(Index1D index_, Index3D shape_) : PeriodIndex3D(toIndex3D(index_, shape_), shape_) {}

    __host__ __device__ inline auto PeriodIndex3D::operator+(const This& other) const -> This {
        assert(shape == other.shape && "[Error]: Inconsistent dimension");
        Index3D result{};
        for (int i = 0; i < 3; ++i)
            result[i] = (index[i] + other.index[i]) % shape[i];
        return PeriodIndex3D(result, shape);
    }

    __host__ __device__ inline auto PeriodIndex3D::toIndex1D() const -> Index1D {
        return (index[0] * shape[1] + index[1]) * shape[2] + index[2];
    }

    __host__ __device__ inline bool PeriodIndex3D::isInReducedK() const {
        const size_t kSpaceDimZ = FFT<Real<>>::rSizeToKSize(shape[2]);
        return index[2] < kSpaceDimZ;
    }

    __host__ __device__ inline Index3D PeriodIndex3D::toReducedK() const {
        Index3D result = index;
        if (!isInReducedK()) {
            result[0] = (shape[0] - result[0]) % shape[0];
            result[1] = (shape[1] - result[1]) % shape[1];
            result[2] = shape[2] - result[2];
        }
        return result;
    }

    __host__ __device__ inline void PeriodIndex3D::normalize() {
        for (int i = 0; i < 3; ++i)
            index[i] %= shape[i];
    }

    __host__ __device__ inline void PeriodIndex3D::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        index.swap(obj.index);
        shape.swap(obj.shape);
    }

    __host__ __device__ inline Index3D PeriodIndex3D::toIndex3D(Index1D index, Index3D shape) {
        Index3D result{};
        result[0] = index / (shape[1] * shape[2]);
        const Index1D temp = index % (shape[1] * shape[2]);
        result[1] = temp / shape[2];
        result[2] = temp % shape[2];
        return result;
    }
}
