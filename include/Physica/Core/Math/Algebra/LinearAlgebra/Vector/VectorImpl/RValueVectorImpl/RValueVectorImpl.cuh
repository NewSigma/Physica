/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Utils/CUDA/PlainStruct.h>
#include <Physica/Core/Parallel/CUDAContext.cuh>

namespace Physica::Core {
    template<class Derived>
    template<class OtherDerived>
    __device__ void device_obj<RValueVector<Derived>>::assignTo(device_obj<LValueVector<OtherDerived>>& target_) const {
        using OtherScalar = typename OtherDerived::ScalarType;
        constexpr size_t OtherSize = Traits<OtherDerived>::SizeAtCompile;
        static_assert(SizeAtCompile == Dynamic || OtherSize == Dynamic || SizeAtCompile == OtherSize, "[Error]: Size mismatch between two vector");

        auto& target = target_.getDerived();
        assert(target.getLength() == getLength() && "[Error]: Size mismatch between two vector");
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = OtherScalar(calc(i));
    }

    template<class Derived>
    __host__ __device__ inline device_obj<TransposeVector<Derived>> device_obj<RValueVector<Derived>>::transpose() const noexcept {
        return device_obj<TransposeVector<Derived>>(*this);
    }

    template<class Derived>
    __device__ inline typename device_obj<RValueVector<Derived>>::RealType device_obj<RValueVector<Derived>>::norm() const {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    __device__ inline typename device_obj<RValueVector<Derived>>::RealType device_obj<RValueVector<Derived>>::squaredNorm() const {
        auto result = RealType(0);
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::max() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::min() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::sum() const {
        assert(getLength() != 0);
        auto result = ScalarType(0);
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i);
        return result;
    }

    template<class Derived>
    template<class OtherDerived>
    __device__ inline device_obj<CrossProduct<Derived, OtherDerived>>
    device_obj<RValueVector<Derived>>::crossProduct(const device_obj<RValueVector<OtherDerived>>& v) const noexcept {
        return device_obj<CrossProduct<Derived, OtherDerived>>(*this, v);
    }
}
