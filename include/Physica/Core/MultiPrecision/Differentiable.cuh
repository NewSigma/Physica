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
    namespace Internal {
        template<class T, DiffMode M>
        class Traits<Core::device_obj<Differentiable<T, M>>> {
            static_assert(!T::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
            static_assert(!Utils::is_device_obj<T>::value, "[Error]: Nested device_obj<> is not allowed");
            using RealT = typename T::RealType;
            using ComplexT = typename T::ComplexType;
        public:
            using ScalarType = device_obj<Differentiable<T, M>>;
            using RealType = device_obj<Differentiable<RealT, M>>;
            using ComplexType = device_obj<Differentiable<ComplexT, M>>;
            using TrivialType = typename T::TrivialType;
            using PlainScalar = T;
            constexpr static ScalarOption Option = T::Option;
            constexpr static bool isComplex = T::isComplex;
            constexpr static bool isDifferentiable = true;
            constexpr static DiffMode Mode = M;
        };
    }

    template<class ScalarType>
    class device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
            : public ScalarBase<device_obj<Differentiable<ScalarType, DiffMode::Reverse>>> {
        using host_obj = Differentiable<ScalarType, DiffMode::Reverse>;
        using This = device_obj<host_obj>;
        using Base = ScalarBase<This>;
        using TracerType = device_obj<DiffTracer<ScalarType>>;

        ScalarType* __restrict pValue;
        ScalarType* __restrict pGrad;
    public:
        device_obj() = default;
        device_obj(double d) : device_obj(ScalarType(d)) {}
        device_obj(ScalarType s);
        __host__ __device__ device_obj(const ScalarType* pValue_, const ScalarType* pGrad_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        [[nodiscard]] __device__ explicit operator float() const { return float(getValue()); }
        [[nodiscard]] __device__ explicit operator double() const { return double(getValue()); }
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const;
        [[nodiscard]] inline device_obj operator-() const;
        /* Operations */
        void toHostAsync_value(ScalarType& value) const;
        void toHostAsync_grad(ScalarType& grad) const;
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
        __device__ void setValue(const ScalarType& x) { *pValue = x; }
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator+(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2);

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator-(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2);

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator*(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2);

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator/(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2);
}

#include "DifferentiableImpl/DifferentiableImpl.cuh"
#include "DifferentiableImpl/ElementaryFunction.cuh"
