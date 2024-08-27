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
    class device_obj<VectorExpr<ExpressionType::Reciprocal, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Reciprocal, VectorType>>> {
        using host_obj = VectorExpr<ExpressionType::Reciprocal, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        const DeviceVector& exp;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(exp_.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return reciprocal(exp.calc(index));
        }
        [[nodiscard]] __device__ size_t getLength() const {
            return exp.getLength();
        }
    };

    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Relu, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Relu, VectorType>>> {
        using host_obj = VectorExpr<ExpressionType::Relu, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        const DeviceVector& exp;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(exp_.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return relu(exp.calc(index));
        }
        [[nodiscard]] __device__ size_t getLength() const {
            return exp.getLength();
        }
    };

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
    [[nodiscard]] __device__ inline auto reciprocal(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Reciprocal, VectorType>>(v);
    }

    template<class VectorType>
    [[nodiscard]] __device__ inline auto relu(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Relu, VectorType>>(v);
    }

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
#include "VectorExprImpl/VectorSquare.cuh"
