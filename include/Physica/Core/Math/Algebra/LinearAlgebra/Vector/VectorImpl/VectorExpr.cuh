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
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class device_obj<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        using DeviceVector = device_obj<VectorType>;
    public:
        using host_obj = VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        Physica::PlainStruct<const DeviceVector> v;
        AnyScalar s;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& v_, AnyScalar s_)
                : v(asStruct(v_.getDerived())), s(s_) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return ScalarType(v.getDerived().calc(index)) / ScalarType(s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getDerived().getLength(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Reciprocal, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Reciprocal, VectorType>>> {
    public:
        using host_obj = VectorExpr<ExpressionType::Reciprocal, VectorType>;
    private:
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        Physica::PlainStruct<const DeviceVector> exp;
    public:
        device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(asStruct(exp_.getDerived())) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return reciprocal(exp.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return exp.getDerived().getLength();
        }
    };

    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Relu, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Relu, VectorType>>> {
    public:
        using host_obj = VectorExpr<ExpressionType::Relu, VectorType>;
    private:
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        Physica::PlainStruct<const DeviceVector> exp;
    public:
        device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(asStruct(exp_.getDerived())) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return relu(exp.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return exp.getDerived().getLength();
        }
    };

    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Exp, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Exp, VectorType>>> {
    public:
        using host_obj = VectorExpr<ExpressionType::Exp, VectorType>;
    private:
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        Physica::PlainStruct<const DeviceVector> v;
    public:
        device_obj(const device_obj<RValueVector<VectorType>>& v_) : v(asStruct(v_.getDerived())) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return exp(v.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return v.getDerived().getLength();
        }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>>
    operator/(const device_obj<RValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) noexcept {
        return {v.getDerived(), s.getDerived()};
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    [[nodiscard]] inline auto reciprocal(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Reciprocal, VectorType>>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline auto relu(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Relu, VectorType>>(v);
    }

    template<class VectorType>
    inline auto exp(const device_obj<RValueVector<VectorType>>& v) noexcept {
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
#include "VectorExprImpl/VectorSquare.cuh"
