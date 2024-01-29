/*
 * Copyright 2024 WeiBo He.
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

#include "DifferentiableImpl/DiffTracer.cuh"
#include "Differentiable.h"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<Differentiable<ScalarType, DiffMode::Reverse>> {
        using host_obj = Differentiable<ScalarType, DiffMode::Reverse>;
        using This = device_obj<host_obj>;

        ScalarType* __restrict pValue;
        ScalarType* __restrict pGrad;
    public:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        [[nodiscard]] __device__ explicit operator float() const { return float(getValue()); }
        [[nodiscard]] __device__ explicit operator double() const { return double(getValue()); }
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const;
        /* Operations */
        __host__ __device__ inline void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ ScalarType* value_ptr() const noexcept { return pValue; }
        [[nodiscard]] __host__ __device__ ScalarType* grad_ptr() const noexcept { return pGrad; }
        [[nodiscard]] __device__ ScalarType& getValue() noexcept { return *pValue; }
        [[nodiscard]] __device__ ScalarType& getGrad() noexcept { return *pGrad; }
        [[nodiscard]] __device__ const ScalarType& getValue() const noexcept { return *pValue; }
        [[nodiscard]] __device__ const ScalarType& getGrad() const noexcept { return *pGrad; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return getValue().isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return getValue().isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return getValue().isNegative(); }
        /* Setters */
        __deivce__ void setValue(const ScalarType& x) { *pValue = x; }
    };

    template<class ScalarType>
    __host__ __device__ inline bool device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (pGrad == other.pGrad) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType>
    __host__ __device__ inline void device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(pGrad, obj.pGrad);
    }
}
