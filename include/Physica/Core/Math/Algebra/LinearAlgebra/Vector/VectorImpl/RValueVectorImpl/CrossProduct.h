/*
 * Copyright 2021-2024 Weibo He.
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

#include "../RValueVector.h"

namespace Physica {
    template<Vector V1, Vector V2>
    class CrossProduct : public RValueVector<CrossProduct<V1, V2>> {
        using This = CrossProduct<V1, V2>;
        using Base = RValueVector<This>;
        using Base::isReverseDiff;
        using typename Base::ScalarType;
    private:
        const V1& v1;
        const V2& v2;
    public:
        CrossProduct(const V1& v1_, const V2& v2_) : v1(v1_), v2(v2_) {
            assert(v1.getLength() == 3);
            assert(v2.getLength() == 3);
        }
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const {
            v[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }

        [[nodiscard]] ScalarType calc(size_t index) const;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return 3; }
    };

    template<Vector V1, Vector V2>
    CrossProduct<V1, V2>::ScalarType CrossProduct<V1, V2>::calc(size_t index) const {
        assert(index < getLength());
        switch (index) {
        case 0:
            return v1[1] * v2[2] - v1[2] * v2[1];
        case 1:
            return v1[2] * v2[0] - v1[0] * v2[2];
        case 2:
            return v1[0] * v2[1] - v1[1] * v2[0];
        default:
            unreachable();
        }
    }
}

namespace Physica {
    template<Vector V1, Vector V2>
    class Traits<CrossProduct<V1, V2>> {
        static_assert((Traits<V1>::SizeAtCompile == 3 || Traits<V1>::SizeAtCompile == Dynamic) &&
                      (Traits<V2>::SizeAtCompile == 3 || Traits<V2>::SizeAtCompile == Dynamic),
                      "CrossProduct can apply on 3-dim vectors only");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = 3;
    };
}
