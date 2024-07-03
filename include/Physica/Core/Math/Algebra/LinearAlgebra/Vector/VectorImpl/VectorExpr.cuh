/*
 * Copyright 2023-2024 WeiBo He.
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
    //////////////////////////////////////Add//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class device_obj<VectorExpr<ExpressionType::Add, VectorType1, VectorType2>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Add, VectorType1, VectorType2>>> {
        using DeviceVector1 = device_obj<VectorType1>;
        using DeviceVector2 = device_obj<VectorType2>;
    public:
        using host_obj = VectorExpr<ExpressionType::Add, VectorType1, VectorType2>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        Physica::PlainStruct<const DeviceVector1> v1;
        Physica::PlainStruct<const DeviceVector2> v2;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType1>>& v1_, const device_obj<RValueVector<VectorType2>>& v2_)
                : v1(asStruct(v1_.getDerived())), v2(asStruct(v2_.getDerived())) {
            assert(v1_.getLength() == v2_.getLength());
        }
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(v1.getDerived().calc(index)) + ScalarType(v2.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getDerived().getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class device_obj<VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        using DeviceVector = device_obj<VectorType>;
    public:
        using host_obj = VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>;
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
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return ScalarType(v.getDerived().calc(index)) - ScalarType(s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getDerived().getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class device_obj<VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>>> {
        using DeviceVector1 = device_obj<VectorType1>;
        using DeviceVector2 = device_obj<VectorType2>;
    public:
        using host_obj = VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        Physica::PlainStruct<const DeviceVector1> v1;
        Physica::PlainStruct<const DeviceVector2> v2;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType1>>& v1_, const device_obj<RValueVector<VectorType2>>& v2_)
                : v1(asStruct(v1_.getDerived())), v2(asStruct(v2_.getDerived())) {
            assert(v1_.getLength() == v2_.getLength());
        }
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(v1.getDerived().calc(index)) - ScalarType(v2.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getDerived().getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class device_obj<VectorExpr<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        using DeviceVector = device_obj<VectorType>;
    public:
        using host_obj = VectorExpr<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>;
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
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return ScalarType(v.getDerived().calc(index)) * ScalarType(s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getDerived().getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class device_obj<VectorExpr<ExpressionType::Mul, VectorType1, VectorType2>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Mul, VectorType1, VectorType2>>> {
        using DeviceVector1 = device_obj<VectorType1>;
        using DeviceVector2 = device_obj<VectorType2>;
    public:
        using host_obj = VectorExpr<ExpressionType::Mul, VectorType1, VectorType2>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        Physica::PlainStruct<const DeviceVector1> v1;
        Physica::PlainStruct<const DeviceVector2> v2;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType1>>& v1_, const device_obj<RValueVector<VectorType2>>& v2_)
                : v1(asStruct(v1_.getDerived())), v2(asStruct(v2_.getDerived())) {
            assert(v1_.getLength() == v2_.getLength());
        }
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return v1.getDerived().calc(index) * v2.getDerived().calc(index); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getDerived().getLength(); }
    };
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
    //////////////////////////////////////Add//////////////////////////////////////
    template<class Derived, class OtherDerived>
    [[nodiscard]] inline __host__ __device__ device_obj<VectorExpr<ExpressionType::Add, Derived, OtherDerived>>
            operator+(const device_obj<RValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) noexcept {
        return {v1.getDerived(), v2.getDerived()};
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<ScalarType>>>
    operator-(const device_obj<RValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) noexcept {
        return {v.getDerived(), s.getDerived()};
    }

    template<class Derived, class OtherDerived>
    [[nodiscard]] inline __host__ __device__ device_obj<VectorExpr<ExpressionType::Sub, Derived, OtherDerived>>
            operator-(const device_obj<RValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) noexcept {
        return {v1.getDerived(), v2.getDerived()};
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>>
    operator*(const device_obj<RValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) noexcept {
        return {v.getDerived(), s.getDerived()};
    }

    template<class ScalarType, class VectorType>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>>
    operator*(const ScalarBase<ScalarType>& s, const device_obj<RValueVector<VectorType>>& v) noexcept {
        return v * s.getDerived();
    }
    
    template<class VectorType1, class VectorType2>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExpressionType::Mul, VectorType1, VectorType2>>
    hadamard(const device_obj<RValueVector<VectorType1>>& v1, const device_obj<RValueVector<VectorType2>>& v2) noexcept {
        return {v1.getDerived(), v2.getDerived()};
    }
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
