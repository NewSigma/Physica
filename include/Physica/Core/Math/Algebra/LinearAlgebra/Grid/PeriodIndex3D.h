/*
 * Copyright 2023 WeiBo He.
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

#include "GridImpl/GridBase.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    /**
     * Provide a algebra for Index1D and Index3D in periodic boundary condition
     */
    class PeriodIndex3D {
        using Index3D = typename GridBase::Index3D;
        using Index1D = typename Index3D::ValueType;

        Index3D index;
        Index3D dim;
    public:
        inline PeriodIndex3D(Index3D index_, Index3D dim_);
        inline PeriodIndex3D(Index1D index_, Index3D dim_);
        PeriodIndex3D(const PeriodIndex3D&) = default;
        PeriodIndex3D(PeriodIndex3D&&) noexcept = default;
        ~PeriodIndex3D() = default;
        /* Operators */
        inline PeriodIndex3D& operator=(PeriodIndex3D obj) noexcept;
        [[nodiscard]] operator Index3D() { return index; }
        [[nodiscard]] inline PeriodIndex3D operator+(const PeriodIndex3D& other) const;
        /* Operations */
        [[nodiscard]] inline Index1D toIndex1D() const;
        [[nodiscard]] inline bool isInReducedK() const;
        [[nodiscard]] inline Index3D toReducedK() const;
        inline void normalize();
        inline void swap(PeriodIndex3D& obj) noexcept;
    private:
        static Index3D toIndex3D(Index1D index, Index3D dim);
    };

    inline PeriodIndex3D::PeriodIndex3D(Index3D index_, Index3D dim_) : index(index_), dim(dim_) {
        for (int i = 0; i < 3; ++i)
            assert(dim[i] > 0 && "[Error]: Zero dim is not allowed");
    }

    inline PeriodIndex3D::PeriodIndex3D(Index1D index_, Index3D dim_) : PeriodIndex3D(toIndex3D(index_, dim_), dim_) {}

    inline PeriodIndex3D& PeriodIndex3D::operator=(PeriodIndex3D obj) noexcept {
        swap(obj);
        return *this;
    }

    inline PeriodIndex3D PeriodIndex3D::operator+(const PeriodIndex3D& other) const {
        assert(dim == other.dim && "[Error]: Inconsistent dimention");
        Index3D result{};
        for (int i = 0; i < 3; ++i)
            result[i] = (index[i] + other.index[i]) % dim[i];
        return PeriodIndex3D(result, dim);
    }

    inline typename PeriodIndex3D::Index1D PeriodIndex3D::toIndex1D() const {
        return (index[0] * dim[1] + index[1]) * dim[2] + index[2];
    }

    inline bool PeriodIndex3D::isInReducedK() const {
        const size_t kSpaceDimZ = FFT<Scalar<>>::rSizeToKSize(dim[2]);
        return index[2] < kSpaceDimZ;
    }

    inline typename PeriodIndex3D::Index3D PeriodIndex3D::toReducedK() const {
        Index3D result = index;
        if (!isInReducedK()) {
            result[0] = (dim[0] - result[0]) % dim[0];
            result[1] = (dim[1] - result[1]) % dim[1];
            result[2] = dim[2] - result[2];
        }
        return result;
    }

    inline void PeriodIndex3D::normalize() {
        for (int i = 0; i < 3; ++i)
            index[i] %= dim[i];
    }

    inline void PeriodIndex3D::swap(PeriodIndex3D& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        index.swap(obj.index);
        dim.swap(obj.dim);
    }

    typename PeriodIndex3D::Index3D PeriodIndex3D::toIndex3D(Index1D index, Index3D dim) {
        Index3D result{};
        result[0] = index / (dim[1] * dim[2]);
        const Index1D temp = index % (dim[1] * dim[2]);
        result[1] = temp / dim[2];
        result[2] = temp % dim[2];
        return result;
    }
}
