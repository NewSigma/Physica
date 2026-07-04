/*
 * Copyright 2023-2026 Weibo He.
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

#include <iosfwd>
#include "Physica/Core/Scalar/Real.h"  // IWYU pragma: export
#include "DiffImpl/DiffCoro.h"

namespace Physica {
    template<class ScalarType> class ScalarPtr;
    /**
     * \class Diff provides auto differential support for scalars
     */
    template<Scalar T, DiffMode Mode, int Order>
    class Diff<T, Mode, Order> : public ScalarBase<Diff<T, Mode, Order>>
                               , public std::conditional<Mode == DiffMode::Forward, CRCoro<Diff<T, Mode, Order>>, Empty>::type {
        using This = Diff<T, Mode, Order>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
        using typename Base::RealType;
        using typename Base::GradType;
        using typename Base::MachineType;
        using Base::isComplex;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    private:
        T v;
        GradType g;
    public:
        Diff() = default;
        __host__ __device__ Diff(MachineType x) : This(T(x)) {}
        __host__ __device__ Diff(T v_);
        __host__ __device__ Diff(T v_, GradType g_);
        Diff(DiffCoro<This>) requires(isReverseDiff()) = delete("[Error]: Copying differential coroutines discards the compute graph");
        template<Scalar U>
        __host__ __device__ explicit(T::Prec < U::Prec) Diff(const U& x);
        Diff(const This&) requires(isForwardDiff()) = default;
        Diff(This&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        __host__ __device__ This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ explicit operator float() const { return float(v); }
        [[nodiscard]] __host__ __device__ explicit operator double() const { return double(v); }
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const;
        /* Operations */
        __host__ __device__ T reverse(GradType grad = 1) const noexcept;
        __host__ __device__ void zero_grad();

        [[nodiscard]] This conjugate() const noexcept;
        [[nodiscard]] CoDiff<This> squaredNorm() const noexcept;

        using Base::random_uniform;
        using Base::random_normal;
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        __host__ __device__ void swap(ScalarRef<This>&& ref) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto real_ptr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto* value_ptr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto* grad_ptr(this auto&&) noexcept;

        using Base::value;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto& grad(this auto&&) noexcept;
        template<int MaskOrder>
        [[nodiscard]] __host__ __device__ auto grad_mask() const noexcept;
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept;
        /* Static members */
        [[nodiscard]] static This fromPhase(RealType phase) noexcept;
        template<RNG R>
        [[nodiscard]] static auto random_uniform();
        template<RNG R>
        [[nodiscard]] static auto random_normal();
        template<RNG R>
        [[nodiscard]] static auto random_any(auto& distribution);
        [[nodiscard]] static const H5::DataType& dtype_hdf5() noexcept;
    };

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const T& obj) requires(T::isDiffable()) {
        return os << obj.value();
    }
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order_>
    class Traits<Diff<T, Mode, Order_>> {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
        static_assert(Order_ > 0, "[Error]: Use plain type instead of 0 order differentiable");
        using RealT = T::RealType;
        using ComplexT = T::ComplexType;
    public:
        constexpr static FloatPrec Prec = T::Prec;
        constexpr static int Order = Order_;
        constexpr static bool isComplex = T::isComplex();
        constexpr static bool isForwardDiff = Mode == DiffMode::Forward;
        constexpr static bool isReverseDiff = Mode == DiffMode::Reverse;
    private:
        using GradType = std::conditional<Order == 1, T, Diff<T, Mode, Order - 1>>::type;
    public:
        using ValueType = T;
        using ScalarType = Diff<T, Mode, Order>;
        using PtrTy = ScalarPtr<ScalarType>;
        using ConstPtrTy = ScalarPtr<const ScalarType>;
        using RefTy = ScalarRef<ScalarType>;
        using ConstRefTy = ScalarRef<const ScalarType>;
        using RealType = Diff<RealT, Mode, Order>;
        using ComplexType = Diff<ComplexT, Mode, Order>;
        using MachineType = T::MachineType;
    };
}

namespace std {
    template<Physica::Scalar T, Physica::DiffMode Mode, int Order>
    struct numeric_limits<Physica::Diff<T, Mode, Order>> : public numeric_limits<T> {};

    template<Physica::Scalar T, Physica::DiffMode Mode, int Order>
    struct formatter<Physica::Diff<T, Mode, Order>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Diff<T, Mode, Order>& obj, auto& ctx) const {
            return formatter<T, char>{}.format(obj.value(), ctx);
        }
    };
}

#include "DiffImpl/DiffImpl.h"
#include "DiffImpl/ScalarPtr.h"
#include "DiffImpl/ScalarRef.h"
#include "DiffImpl/MathForward.h"
#include "DiffImpl/MathReverse.h"
#include "DiffImpl/SIMD.h"
