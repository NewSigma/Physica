/*
 * Copyright 2021-2026 Weibo He.
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
    class Cross : public RValueVector<Cross<V1, V2>> {
        using This = Cross<V1, V2>;
        using Base = RValueVector<This>;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<V1> v1;
        decay_rvalue_t<V2> v2;
    public:
        Cross(V1 v1_, V2 v2_);
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;

        [[nodiscard]] T calc(size_t index) const;
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() noexcept { return 3; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return 3; }
    };

    template<Vector V1, Vector V2>
    Cross<V1, V2>::Cross(V1 v1_, V2 v2_) : v1(std::forward<V1>(v1_)), v2(std::forward<V2>(v2_)) {
        constexpr size_t Size1 = std::remove_cvref_t<V1>::getSizeAtCompile();
        constexpr size_t Size2 = std::remove_cvref_t<V2>::getSizeAtCompile();
        static_assert((Size1 == 3 || Size1 == Dynamic) && (Size2 == 3 || Size2 == Dynamic), "[Error]: Cross can apply on 3-dim vectors only");
        assert(v1.getLength() == 3);
        assert(v2.getLength() == 3);
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy>
    void Cross<V1, V2>::assign(Vector auto&& v) const noexcept {
        v[0] = fma(v1[1], v2[2], -v1[2] * v2[1]);
        v[1] = fma(v1[2], v2[0], -v1[0] * v2[2]);
        v[2] = fma(v1[0], v2[1], -v1[1] * v2[0]);
    }

    template<Vector V1, Vector V2>
    auto Cross<V1, V2>::calc(size_t index) const -> T {
        assert(index < getLength());
        switch (index) {
        case 0:
            return fma(v1[1], v2[2], -v1[2] * v2[1]);
        case 1:
            return fma(v1[2], v2[0], -v1[0] * v2[2]);
        case 2:
            return fma(v1[0], v2[1], -v1[1] * v2[0]);
        default:
            unreachable();
        }
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto cross(V1&& v1, V2&& v2) noexcept {
        return Cross<V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

namespace Physica {
    template<Vector V1, Vector V2>
    class Traits<Cross<V1, V2>> {
        using U1 = std::remove_cvref_t<V1>;
        using U2 = std::remove_cvref_t<V2>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename U1::ScalarType, typename U2::ScalarType>::Type;
    };
}
