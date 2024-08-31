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

namespace Physica::Core {
    template<ExpressionType ExprType, class VectorType>
    class device_obj<UnitaryVectorExpr<ExprType, VectorType>> : public device_obj<RValueVector<VectorExpr<ExprType, VectorType>>> {
        using host_obj = UnitaryVectorExpr<ExprType, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ExprType, VectorType>>>;
        using DeviceVector = device_obj<VectorType>;
    public:
        using typename Base::ScalarType;
    private:
        union {
            PlainStruct<const DeviceVector> value;
            const DeviceVector* ptr;
        } expr;
    public:
        __host__ __device__ inline device_obj(const device_obj<RValueVector<VectorType>>& expr_) {
            if constexpr (IsHost())
                expr.value = asStruct(expr_.getDerived());
            else
                expr.ptr = &expr_.getDerived();
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return getExpr<Owner>().template getLength<Owner>();
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceVector& getExpr() const {
            if constexpr (IsHost() || Owner == Side::Host)
                return expr.value.getDerived();
            else
                return *expr.ptr;
        }
    };

    template<ExpressionType ExprType, class VectorType, class T>
    class device_obj<BinaryVectorExpr<ExprType, VectorType, T>>
            : public device_obj<RValueVector<VectorExpr<ExprType, VectorType, T>>> {
        using host_obj = BinaryVectorExpr<ExprType, VectorType, T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ExprType, VectorType, T>>>;
        using DeviceVector = device_obj<VectorType>;
        using DeviceType = typename std::conditional<is_scalar<T>::value, typename T::ScalarType, device_obj<T>>::type;
    public:
        using typename Base::ScalarType;
    private:
        union {
            PlainStruct<const DeviceVector> value;
            const DeviceVector* ptr;
        } lhs;

        union {
            PlainStruct<const DeviceType> value;
            const DeviceType* ptr;
        } rhs;
    public:
        __host__ __device__ inline device_obj(
                const device_obj<RValueVector<VectorType>>& lhs_, const DeviceType& rhs_) {
            if constexpr (!is_scalar<T>::value)
                assert(lhs_.getLength() == rhs_.getLength());
            if constexpr (IsHost()) {
                lhs.value = asStruct(lhs_.getDerived());
                rhs.value = asStruct(rhs_);
            }
            else {
                lhs.ptr = &lhs_.getDerived();
                rhs.ptr = &rhs_;
            }
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getLHS<Owner>().template getLength<Owner>(); }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceVector& getLHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return lhs.value.getDerived();
            else
                return *lhs.ptr;
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceType& getRHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return rhs.value.getDerived();
            else
                return *rhs.ptr;
        }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class device_obj<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        using host_obj = VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>;
        using Base = device_obj<RValueVector<host_obj>>;
        using This = device_obj<host_obj>;
        using DeviceVector = device_obj<VectorType>;
    public:
        using typename Base::ScalarType;
    private:
        const DeviceVector& v;
        AnyScalar s;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType>>& v_, AnyScalar s_) : v(v_.getDerived()), s(s_) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return ScalarType(v.calc(index)) / ScalarType(s); }
        [[nodiscard]] __device__ size_t getLength() const { return v.getLength(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Exp, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Exp, VectorType>>> {
        using host_obj = VectorExpr<ExpressionType::Exp, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        const DeviceVector& v;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType>>& v_) : v(v_.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return exp(v.calc(index));
        }
        [[nodiscard]] __device__ size_t getLength() const {
            return v.getLength();
        }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] __device__ inline device_obj<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>>
    operator/(const device_obj<RValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) noexcept {
        return {v.getDerived(), s.getDerived()};
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    [[nodiscard]] __device__ inline auto exp(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Exp, VectorType>>(v);
    }
}

namespace Physica {
    using namespace Core;

    template<ExpressionType Type, class Exp1, class Exp2>
    class Traits<Core::device_obj<VectorExpr<Type, Exp1, Exp2>>> : public Traits<VectorExpr<Type, Exp1, Exp2>> {};
}

#include "VectorExprImpl/VectorAdd.cuh"
#include "VectorExprImpl/VectorSub.cuh"
#include "VectorExprImpl/VectorMul.cuh"
#include "VectorExprImpl/VectorMinus.cuh"
#include "VectorExprImpl/Reciprocal.cuh"
#include "VectorExprImpl/Relu.cuh"
#include "VectorExprImpl/VectorSquare.cuh"
