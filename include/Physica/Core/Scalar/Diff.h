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

#include <iosfwd>
#include "Physica/Core/Scalar/Real.h"  // IWYU pragma: export
#include "Physica/Core/Scalar/ExprType.h"
#include "DiffImpl/DiffCoro.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    class Diff<T, Mode, Order> : public ScalarBase<Diff<T, Mode, Order>>
                               , public std::conditional<Mode == DiffMode::Forward, CRCoro<Diff<T, Mode, Order>>, PlainStruct<void>>::type {
        using This = Diff<T, Mode, Order>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
        using typename Base::GradType;
        using typename Base::MachineType;
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
        Diff(DiffCoro<This>) requires(isReverseDiff) = delete;
        template<Scalar U>
        __host__ __device__ explicit(T::Option < U::Option) Diff(const U& x) requires(!ReverseDiff<U>);
        Diff(const This&) requires(isForwardDiff) = default;
        Diff(This&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        __host__ __device__ This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ explicit operator float() const { return float(v); }
        [[nodiscard]] __host__ __device__ explicit operator double() const { return double(v); }
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const;
        /* Operations */
        template<int MaskOrder>
        __host__ __device__ auto mask() const noexcept;

        __host__ __device__ T reverse(GradType grad_ = 1) const noexcept;
        __host__ __device__ inline void zero_grad();

        [[nodiscard]] auto conjugate() const;
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        __host__ __device__ void swap(ScalarRef<This>&& ref) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T* value_ptr() noexcept { return &v; }
        [[nodiscard]] __host__ __device__ GradType* grad_ptr() noexcept { return &g; }

        using Base::value;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ inline auto& grad() noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ inline const auto& grad() const noexcept;
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] inline static auto random_uniform();
        template<RNG R>
        [[nodiscard]] inline static auto random_normal();
        template<RNG R, class Distribution>
        [[nodiscard]] inline static auto random_any(Distribution& dist);
        [[nodiscard]] static const H5::DataType& getH5DataType();
    };

    template<Scalar T>
    inline std::ostream& operator<<(std::ostream& os, const T& obj) requires(T::isDiffable) {
        return os << obj.value();
    }
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order_>
    class Traits<Diff<T, Mode, Order_>> {
        static_assert(!T::isDiffable, "[Error]: Nested Diff<> is not allowed");
        static_assert(Order_ > 0, "[Error]: Use plain type instead of 0 order differentiable");
        using RealT = T::RealType;
        using ComplexT = T::ComplexType;
    public:
        constexpr static ScalarOption Option = T::Option;
        constexpr static int Order = Order_;
        constexpr static bool isComplex = T::isComplex;
        constexpr static bool isForwardDiff = Mode == DiffMode::Forward;
        constexpr static bool isReverseDiff = Mode == DiffMode::Reverse;
    private:
        using GradType = std::conditional<Order == 1, T, Diff<T, Mode, Order - 1>>::type;
    public:
        using ValueType = T;
        using ScalarType = Diff<T, Mode, Order>;
        using PtrTy = ScalarPtr<ScalarType>;
        using ConstPtrTy = const ScalarPtr<ScalarType>;
        using RefTy = ScalarRef<ScalarType>;
        using ConstRefTy = const ScalarRef<ScalarType>;
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

        auto format(const Physica::Diff<T, Mode, Order>& obj, std::format_context& ctx) const {
            return formatter<T, char>{}.format(obj, ctx);
        }
    };
}

#include "DiffImpl/DiffImpl.h"
#include "DiffImpl/ScalarPtr.h"
#include "DiffImpl/ScalarRef.h"
#include "DiffImpl/MathForward.h"
#include "DiffImpl/MathReverse.h"
#include "DiffImpl/SIMD.h"
