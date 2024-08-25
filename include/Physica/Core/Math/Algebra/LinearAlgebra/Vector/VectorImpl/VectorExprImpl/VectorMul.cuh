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

namespace Physica::Core {
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
}
