/*
 * Copyright 2024-2025 Weibo He.
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

#include "../VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprType::Exp, V> : public UnitaryVectorExpr<ExprType::Exp, V> {
        using This = VectorExpr<ExprType::Exp, V>;
        using Base = UnitaryVectorExpr<ExprType::Exp, V>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        void reverse(const auto& grad) const noexcept;
    };

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Exp, V>::assign(Vector auto&& v) const noexcept {
        using V1 = std::remove_cvref_t<decltype(v)>;
        constexpr size_t Length = std::max(Base::SizeAtCompile, V1::SizeAtCompile);
        constexpr bool SmallVector = 0 < Length && Length <= 32;
        if constexpr (Internal::EnableMKL<V, V1>::value && !SmallVector && T::Prec == Float64) {
            if (Base::getLength() <= 32)
                Base::template assign_base<P>(v);
            else
                assign_mkl(v);
        }
        else
            Base::template assign_base<P>(v);
    }

    template<Vector V>
    auto VectorExpr<ExprType::Exp, V>::calc(size_t index) const -> CoDiff<T> {
        return exp(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprType::Exp, V>::calc_value(size_t index) const -> Tv {
        return exp(Base::getExpr().calc_value(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Exp, V>::packet(size_t index) const {
        return exp(Base::getExpr().template packet<Pack>(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Exp, V>::packetPartial(size_t index, size_t count) const {
        return exp(Base::getExpr().template packetPartial<Pack>(index, count)).cutoff(count);
    }

    template<Vector V>
    void VectorExpr<ExprType::Exp, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        using U = decltype(grad);
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<U>)
            expr.reverse(exp(expr.values()) * grad.value());
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            expr.reverse(hadamard(exp(expr.values()), grad.values()));
        }
    }

    template<Vector V>
    [[nodiscard]] auto exp(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Exp, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Exp.h"
#endif
