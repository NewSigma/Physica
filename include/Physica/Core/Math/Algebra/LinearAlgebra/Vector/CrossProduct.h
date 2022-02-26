/*
 * Copyright 2021 WeiBo He.
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

#include "Vector.h"

namespace Physica::Core {
    template<class AnyVector1, class AnyVector2> class CrossProduct;

    namespace Internal {
        template<class AnyVector1, class AnyVector2>
        class Traits<CrossProduct<AnyVector1, AnyVector2>> {
        public:
            using ScalarType = typename Internal::BinaryScalarOpReturnType<typename AnyVector1::ScalarType,
                                                                           typename AnyVector2::ScalarType>::Type;
            constexpr static size_t SizeAtCompile = 3;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class AnyVector1, class AnyVector2>
    class CrossProduct : public RValueVector<CrossProduct<AnyVector1, AnyVector2>> {
        static_assert((Internal::Traits<AnyVector1>::SizeAtCompile == 3 || Internal::Traits<AnyVector1>::SizeAtCompile == Dynamic) &&
                      (Internal::Traits<AnyVector2>::SizeAtCompile == 3 || Internal::Traits<AnyVector2>::SizeAtCompile == Dynamic),
                      "CrossProduct can apply on 3-dim vectors only");
        using ReturnType = Vector<typename Internal::Traits<CrossProduct>::ScalarType, 3, 3>;
    public:
        const AnyVector1& v1;
        const AnyVector2& v2;
    public:
        CrossProduct(const RValueVector<AnyVector1>& v1_, const RValueVector<AnyVector2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == 3);
            assert(v2.getLength() == 3);
        }
        /* Operations */
        [[nodiscard]] ReturnType compute() const noexcept {
            ReturnType result{};
            result[0] = v1[1] * v2[2] - v1[2] * v2[1];
            result[1] = v1[2] * v2[0] - v1[0] * v2[2];
            result[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return result;
        }

        template<class OtherDerived>
        void assignTo(LValueVector<OtherDerived>& v) const {
            v[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return 3; }
    };
}
