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
    template<Vector T, Scalar U>
    class VectorExpr<ExprType::Mul, T, U>
            : public BinaryVectorExpr<ExprType::Mul, T, U> {
        using Base = BinaryVectorExpr<ExprType::Mul, T, U>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& v) const;

        [[nodiscard]] ScalarType calc(size_t index) const { return ScalarType(Base::getLHS().calc(index)) * ScalarType(Base::getRHS()); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) * Base::getRHS();
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) * Base::getRHS();
        }

        [[nodiscard]] ScalarType sum() const { return Base::getLHS().sum() * Base::getRHS(); }
    };

    template<Vector T, Scalar U>
    template<LVector V, class Executor>
    inline void VectorExpr<ExprType::Mul, T, U>::assignTo(V& v) const {
        constexpr bool FastAssign = Traits<T>::FastAssign;
        if constexpr (FastAssign) {
            Base::getLHS().template assignTo<V, Executor>(v);
            v *= Base::getRHS();
        }
        else
            Base::assignTo(v);
    }

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Mul, T1, T2>
            : public BinaryVectorExpr<ExprType::Mul, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Mul, T1, T2>;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) * Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) * Base::getRHS().template packetPartial<AnyPacket>(index, count);
        }
    };

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator*(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::Mul, T, U>(v, x);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator*(const U& x, const T& v) noexcept {
        return v * x;
    }
    
    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto hadamard(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::Mul, T1, T2>(v1, v2);
    }
}
