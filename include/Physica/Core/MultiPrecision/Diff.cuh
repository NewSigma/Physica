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

#include "DiffImpl/DiffTracer.cuh"
#include "Diff.h"

namespace Physica::Core {
    template<Scalar T, int Order>
    class device_obj<Diff<T, DiffMode::Reverse, Order>>
            : public ScalarBase<device_obj<Diff<T, DiffMode::Reverse, Order>>> {
        static_assert(Order == 1, "[Error]: High order autodiff is not implemented");
        using host_obj = Diff<T, DiffMode::Reverse, Order>;
        using This = device_obj<host_obj>;
        using Base = ScalarBase<This>;
        using ValuePtr = T* __restrict;
    public:
        using TracerType = device_obj<typename host_obj::TracerType>;
        using typename Base::GradType;
    private:
        ValuePtr pValue;
        GradType pGrad;
    public:
        device_obj() = default;
        device_obj(double d) : device_obj(T(d)) {}
        device_obj(T s);
        __host__ __device__ device_obj(const T* pValue_, const T* pGrad_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        __host__ __device__ inline This& operator=(std::nullptr_t) noexcept;
        [[nodiscard]] __device__ explicit operator float() const { return float(getValue()); }
        [[nodiscard]] __device__ explicit operator double() const { return double(getValue()); }
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const;
        [[nodiscard]] __host__ __device__ inline bool operator==(std::nullptr_t) const noexcept;
        [[nodiscard]] __host__ __device__ inline bool operator!=(std::nullptr_t) const noexcept { return !((*this) == nullptr); }
        [[nodiscard]] inline This operator-() const;
        /* Operations */
        void reverse();
        void reverse_to(This to);
        [[nodiscard]] T toHost_value() const;
        [[nodiscard]] T toHost_grad() const;
        void toHostAsync_value(T& value) const;
        void toHostAsync_grad(T& grad) const;
        __host__ __device__ inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T* value_ptr() const noexcept { return pValue; }
        [[nodiscard]] __host__ __device__ T* grad_ptr() const noexcept { return pGrad; }
        [[nodiscard]] __device__ T& getValue() noexcept { return *pValue; }
        [[nodiscard]] __device__ T& getGrad() noexcept { return *pGrad; }
        [[nodiscard]] __device__ const T& getValue() const noexcept { return *pValue; }
        [[nodiscard]] __device__ const T& getGrad() const noexcept { return *pGrad; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return getValue().isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return getValue().isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return getValue().isNegative(); }
        /* Setters */
        __device__ void setValue(T value) { *pValue = value; }
        __device__ void setGrad(T grad) { *pGrad = grad; }
    };
    ////////////////////////////////////////////////////////////
    template<Scalar T, int Order>
    [[nodiscard]] inline device_obj<Diff<T, DiffMode::Reverse, Order>>
    operator+(const device_obj<Diff<T, DiffMode::Reverse, Order>>& s1,
              const device_obj<Diff<T, DiffMode::Reverse, Order>>& s2);

    template<Scalar T, int Order>
    [[nodiscard]] inline device_obj<Diff<T, DiffMode::Reverse, Order>>
    operator-(const device_obj<Diff<T, DiffMode::Reverse, Order>>& s1,
              const device_obj<Diff<T, DiffMode::Reverse, Order>>& s2);

    template<Scalar T, int Order>
    [[nodiscard]] inline device_obj<Diff<T, DiffMode::Reverse, Order>>
    operator*(const device_obj<Diff<T, DiffMode::Reverse, Order>>& s1,
              const device_obj<Diff<T, DiffMode::Reverse, Order>>& s2);

    template<Scalar T, int Order>
    [[nodiscard]] inline device_obj<Diff<T, DiffMode::Reverse, Order>>
    operator/(const device_obj<Diff<T, DiffMode::Reverse, Order>>& s1,
              const device_obj<Diff<T, DiffMode::Reverse, Order>>& s2);
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order_>
    class Traits<Core::device_obj<Diff<T, Mode, Order_>>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj<> is not allowed");
        using RealT = typename T::RealType;
        using ComplexT = typename T::ComplexType;
    public:
        constexpr static ScalarOption Option = T::Option;
        constexpr static int Order = Order_;
        constexpr static bool isComplex = T::isComplex;
        constexpr static bool isDifferentiable = true;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = true;

        using ValueType = T;
        using ScalarType = Core::device_obj<Diff<T, Mode, Order>>;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = Core::device_obj<Diff<RealT, Mode, Order>>;
        using ComplexType = Core::device_obj<Diff<ComplexT, Mode, Order>>;
        using MachineType = typename T::MachineType;
    };
}

#include "DiffImpl/DiffImpl.cuh"
#include "DiffImpl/Math.cuh"
