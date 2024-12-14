/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/MultiPrecision/Real.h"  // IWYU pragma: export
#include "Physica/Core/MultiPrecision/ExprType.h"
#include "DiffImpl/CoDiff.h"

namespace Physica::Core {
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
        T value;
        GradType grad;
    public:
        Diff() = default;
        Diff(MachineType x) : This(T(x)) {}
        Diff(T value_);
        Diff(T value_, GradType grad_);
        template<Scalar U>
        explicit Diff(const U& x);
        Diff(const This&) = default;
        Diff(This&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] explicit operator float() const { return float(value); }
        [[nodiscard]] explicit operator double() const { return double(value); }
        [[nodiscard]] inline bool operator==(const This& other) const;
        /* Operations */
        template<int MaskOrder>
        auto& mask() noexcept;
        template<int MaskOrder>
        const auto& mask() const noexcept;

        T reverse(GradType grad_ = 1) const noexcept;
        inline void zero_grad();

        [[nodiscard]] auto conjugate() const;
        void swap(This& __restrict obj) noexcept;
        void swap(ScalarRef<This>&& ref) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T* value_ptr() noexcept { return &value; }
        [[nodiscard]] __host__ __device__ GradType* grad_ptr() noexcept { return &grad; }
        [[nodiscard]] T& getValue() noexcept { return value; }
        [[nodiscard]] const T& getValue() const noexcept { return value; }
        template<int GradOrder = 1>
        [[nodiscard]] inline Base::template GradRtnTy<GradOrder>& getGrad() noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] inline const Base::template GradRtnTy<GradOrder>& getGrad() const noexcept;
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] inline static auto random_uniform();
        template<RandomGenerator R>
        [[nodiscard]] inline static auto random_normal();
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static auto random_any(Distribution& dist);
        [[nodiscard]] static const H5::DataType& getH5DataType();
    };

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator+(T&& x, U&& y) requires(Diffable<T>);

    template<Scalar T, Scalar U> 
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator-(T&& x, U&& y) requires(Diffable<T>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator*(T&& x, U&& y) requires(Diffable<T>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator/(T&& x, U&& y) requires(Diffable<T>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator+(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator-(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator*(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>);

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator/(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>);

    template<Scalar T>
    [[nodiscard]] inline CoDiff<T> operator-(T&& x) requires(Diffable<T>);

    template<Scalar T>
    inline std::ostream& operator<<(std::ostream& os, const T& obj) requires(T::isDiffable) {
        return os << obj.getValue();
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
        constexpr static ExprType Source = ExprType::Set;
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
    template<Physica::Core::Scalar T, Physica::Core::DiffMode Mode, int Order>
    struct numeric_limits<Physica::Core::Diff<T, Mode, Order>> : public numeric_limits<T> {};
}

#include "DiffImpl/DiffImpl.h"
#include "DiffImpl/ScalarPtr.h"
#include "DiffImpl/ScalarRef.h"
#include "DiffImpl/Math.h"
#include "DiffImpl/SIMD.h"
