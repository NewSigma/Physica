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
    class device_obj<VectorExpr<ExprType::Add, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<BinaryVectorExpr<ExprType::Add, VectorType, ScalarBase<AnyScalar>>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Add, VectorType, ScalarBase<AnyScalar>>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(Base::template getLHS<Owner>().template calc<Owner>(index))
                 + ScalarType(Base::template getRHS<Owner>());
        }

        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ AnyPacket packet(size_t index) const {
            return Base::template getLHS<Owner>().template packet<AnyPacket, Owner>(index)
                 + AnyPacket(Base::template getRHS<Owner>());
        }

        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::template getLHS().template packetPartial<AnyPacket, Owner>(index, count)
                 + AnyPacket(Base::template getRHS());
        }
    };

    template<class VectorType1, class VectorType2>
    class device_obj<VectorExpr<ExprType::Add, VectorType1, VectorType2>>
            : public device_obj<BinaryVectorExpr<ExprType::Add, VectorType1, VectorType2>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Add, VectorType1, VectorType2>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(Base::template getLHS<Owner>().template calc<Owner>(index))
                 + ScalarType(Base::template getRHS<Owner>().template calc<Owner>(index));
        }

        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ AnyPacket packet(size_t index) const {
            return Base::template getLHS<Owner>().template packet<AnyPacket, Owner>(index)
                 + Base::template getRHS<Owner>().template packet<AnyPacket, Owner>(index);
        }

        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::template getLHS<Owner>().template packetPartial<AnyPacket, Owner>(index, count)
                 + Base::template getRHS<Owner>().template packetPartial<AnyPacket, Owner>(index, count);
        }
    };

    template<class VectorType, class AnyScalar>
    [[nodiscard]] __host__ __device__ inline auto operator+(
            const device_obj<RValueVector<VectorType>>& v, const ScalarBase<AnyScalar>& s) noexcept {
        return device_obj<VectorExpr<ExprType::Add, VectorType, ScalarBase<AnyScalar>>>(v.getDerived(), s.getDerived());
    }

    template<class Derived, class OtherDerived>
    [[nodiscard]] __host__ __device__ inline auto operator+(
            const device_obj<RValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) noexcept {
        return device_obj<VectorExpr<ExprType::Add, Derived, OtherDerived>>(v1.getDerived(), v2.getDerived());
    }
}
