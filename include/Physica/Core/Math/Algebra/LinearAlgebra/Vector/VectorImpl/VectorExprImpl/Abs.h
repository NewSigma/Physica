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
    class VectorExpr<ExprType::Abs, V> : public UnitaryVectorExpr<ExprType::Abs, V> {
        using This = VectorExpr<ExprType::Abs, V>;
        using Base = UnitaryVectorExpr<ExprType::Abs, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        constexpr static bool isComplexV = std::remove_cvref_t<V>::isComplex;
        constexpr static bool isReverseDiff = Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return abs(Base::getExpr().calc(index)); }

        [[nodiscard]] Tv calc_value(size_t index) const { return abs(Base::getExpr().calc_value(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            if constexpr (isComplexV)
                return sqrt(Base::getExpr().squaredNorms().template packet<Pack>(index));
            else
                return abs(Base::getExpr().template packet<Pack>(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            if constexpr (isComplexV)
                return sqrt(Base::getExpr().squaredNorms().template packetPartial<Pack>(index, count));
            else
                return abs(Base::getExpr().template packetPartial<Pack>(index, count));
        }

        void reverse(const auto& grad) const noexcept requires(isReverseDiff);

        T max() const {
            if constexpr (isComplexV)
                return sqrt(Base::getExpr().squaredNorms().max());
            else
                return Base::max();
        }
    };

    template<Vector V>
    void VectorExpr<ExprType::Abs, V>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        using U = decltype(grad);
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<U>)
            expr.reverse(unit(expr.values()) * grad);
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            expr.reverse(hadamard(unit(expr.values()), grad));
        }
    }

    template<Vector V>
    [[nodiscard]] inline auto abs(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Abs, V&&>(std::forward<V>(v));
    }
}
