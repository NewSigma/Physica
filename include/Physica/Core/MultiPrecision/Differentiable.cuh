/*
 * Copyright 2024 Weibo He.
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
    template<class ScalarType, unsigned int Order>
    class device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
            : public ScalarBase<device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>> {
        static_assert(Order == 1, "[Error]: High order autodiff is not implemented");
        using host_obj = Differentiable<ScalarType, DiffMode::Reverse, Order>;
        using This = device_obj<host_obj>;
        using Base = ScalarBase<This>;
        using ReducedType = device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order - 1>>;
        using ValueType = ScalarType* __restrict;
        using GradType = typename std::conditional<Order == 1, ValueType, ReducedType>::type;
    public:
        using TracerType = device_obj<typename host_obj::TracerType>;
    private:
        ValueType pValue;
        GradType pGrad;
    public:
        device_obj() = default;
        device_obj(double d) : device_obj(ScalarType(d)) {}
        device_obj(ScalarType s);
        __host__ __device__ device_obj(const ScalarType* pValue_, const ScalarType* pGrad_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] __device__ explicit operator float() const { return float(getValue()); }
        [[nodiscard]] __device__ explicit operator double() const { return double(getValue()); }
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const;
        [[nodiscard]] inline This operator-() const;
        /* Operations */
        void reverse();
        void reverse_to(This to);
        [[nodiscard]] ScalarType toHost_value() const;
        [[nodiscard]] ScalarType toHost_grad() const;
        void toHostAsync_value(ScalarType& value) const;
        void toHostAsync_grad(ScalarType& grad) const;
        __host__ __device__ inline void swap(This& __restrict obj) noexcept;
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
        __device__ void setValue(ScalarType value) { *pValue = value; }
        __device__ void setGrad(ScalarType grad) { *pGrad = grad; }
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator+(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2);

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator-(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2);

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator*(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2);

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator/(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2);
}

namespace Physica {
    using namespace Core;

    template<class T, DiffMode M, unsigned int Order_>
    class Traits<Core::device_obj<Differentiable<T, M, Order_>>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        static_assert(!Utils::is_device_obj<T>::value, "[Error]: Nested device_obj<> is not allowed");
        using RealT = typename T::RealType;
        using ComplexT = typename T::ComplexType;
    public:
        using PlainScalar = T;
        constexpr static DiffMode Mode = M;
        constexpr static unsigned int Order = Order_;
        using ScalarType = Core::device_obj<Differentiable<T, M, Order>>;
        using RealType = Core::device_obj<Differentiable<RealT, M, Order>>;
        using ComplexType = Core::device_obj<Differentiable<ComplexT, M, Order>>;
        using TrivialType = typename T::TrivialType;
        constexpr static ScalarOption Option = T::Option;
        constexpr static bool isComplex = T::isComplex;
        constexpr static bool isDifferentiable = true;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = true;
    };
}

#include "DifferentiableImpl/DifferentiableImpl.cuh"
#include "DifferentiableImpl/ElementaryFunction.cuh"
