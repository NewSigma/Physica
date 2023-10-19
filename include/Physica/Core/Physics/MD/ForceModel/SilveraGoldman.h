/*
 * Copyright 2023 WeiBo He.
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

#include "PairModel.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType> class SilveraGoldman;

    namespace Internal {
        template<class T, class U>
        class Traits<SilveraGoldman<T, U>> {
        public:
            using ScalarType = T;
            using PosScalarType = U;
            constexpr static bool IsPotDependOnAtomIndex = false;
        };
    }
    /**
     * Potential that suits para-hydrogen
     * 
     * Reference:
     * [1] I. F. Silvera and V. V. Goldman, J. Chem. Phys. 69, 4209 (1978).
     */
    template<class ScalarType, class PosScalarType>
    class SilveraGoldman final : public PairModel<SilveraGoldman<ScalarType, PosScalarType>> {
        constexpr static double alpha = 1.713;
        constexpr static double beta = 1.5671;
        constexpr static double gamma = 0.00993;
        constexpr static double cutoff = 8.32;
        constexpr static double c6 = 12.14;
        constexpr static double c8 = 215.2;
        constexpr static double c9 = 143.1;
        constexpr static double c10 = 4813.9;
        constexpr static double squaredCutoff = cutoff * cutoff;

        using This = SilveraGoldman<ScalarType, PosScalarType>;
        using Base = PairModel<This>;
        using typename Base::PlainScalar;
    public:
        SilveraGoldman(PlainScalar cutoff_);
        SilveraGoldman(const SilveraGoldman&) = default;
        SilveraGoldman(SilveraGoldman&&) noexcept = default;
        ~SilveraGoldman() = default;
        /* Operators */
        SilveraGoldman& operator=(SilveraGoldman obj) noexcept;
        /* Operations */
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType forceConst_functor(ScalarType r, ScalarType r2) const;
        inline void swap(SilveraGoldman& obj) noexcept;
    };

    template<class ScalarType, class PosScalarType>
    SilveraGoldman<ScalarType, PosScalarType>::SilveraGoldman(PlainScalar cutoff_) : Base(std::move(cutoff_)) {}

    template<class ScalarType, class PosScalarType>
    SilveraGoldman<ScalarType, PosScalarType>& SilveraGoldman<ScalarType, PosScalarType>::operator=(SilveraGoldman obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    inline void SilveraGoldman<ScalarType, PosScalarType>::swap(SilveraGoldman& obj) noexcept {
        Base::swap(obj);
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType SilveraGoldman<ScalarType, PosScalarType>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, ScalarType r2) const {
        ScalarType result = exp(-r2 * PlainScalar(gamma) - r * PlainScalar(beta) + PlainScalar(alpha));
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r8 = square(rep_r4);
        const ScalarType rep_r9 = rep_r8 * rep_r;
        const ScalarType rep_r10 = rep_r6 * rep_r4;
        const ScalarType g = rep_r6 * PlainScalar(c6) + rep_r8 * PlainScalar(c8) - rep_r9 * PlainScalar(c9) + rep_r10 * PlainScalar(c10);

        if (r < PlainScalar(cutoff)) {
            const ScalarType f_cutoff = exp(-square(rep_r * PlainScalar(cutoff) - PlainScalar(1)));
            result -= g * f_cutoff;
        }
        else
            result -= g;
        return result;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType SilveraGoldman<ScalarType, PosScalarType>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, ScalarType r2) const {
        using PacketType = typename Internal::BestPacket<PlainScalar, 4>::Type;
        constexpr bool isSIMDEnabled = !std::is_same<PlainScalar, PacketType>::value;

        const ScalarType factor = r * PlainScalar(gamma * 2) + PlainScalar(beta);
        ScalarType result = exp(-r2 * PlainScalar(gamma) - (r * PlainScalar(beta) - PlainScalar(alpha))) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;
        const ScalarType rep_r6 = rep_r5 * rep_r;

        ScalarType d_g;
        if constexpr (isSIMDEnabled) {
            const auto rep = PacketType(double(rep_r), double(rep_r3), double(rep_r4), double(rep_r5));
            const PacketType rep1 = rep * rep_r6.getValue();
            const PacketType const1(-6 * c6, -8 * c8, 9 * c9, -10 * c10);
            d_g = (rep1 * const1).horizontal_add();
        }
        else {
            d_g = rep_r * PlainScalar(-6 * c6) + rep_r3 * PlainScalar(-8 * c8) + rep_r4 * PlainScalar(9 * c9) + rep_r5 * PlainScalar(-10 * c10);
            d_g *= rep_r6;
        }

        if (r < PlainScalar(cutoff)) {
            const ScalarType g = (rep_r2 * PlainScalar(c6) + rep_r4 * PlainScalar(c8)) * rep_r4 - (rep_r5 * PlainScalar(c9) - rep_r6 * PlainScalar(c10)) * rep_r4;
            const ScalarType f_cutoff = exp(-square(rep_r * PlainScalar(cutoff) - PlainScalar(1)));
            result += (d_g + g * (rep_r3 * PlainScalar(2 * squaredCutoff) - rep_r2 * PlainScalar(2 * cutoff))) * f_cutoff;
        }
        else
            result += d_g;
        return result;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType SilveraGoldman<ScalarType, PosScalarType>::forceConst_functor(ScalarType r, ScalarType r2) const {
        using PacketType = typename Internal::BestPacket<ScalarType, 4>::Type;
        static_assert(!std::is_same<ScalarType, PacketType>::value, "[Error]: SIMD is inavailable, implementation must be revised");

        const ScalarType factor = square(r * (gamma * 2) + beta) - ScalarType(2 * gamma);
        ScalarType result = exp(-r2 * gamma - (r * beta - alpha)) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;
        const PacketType rep(rep_r.getTrivial(), rep_r3.getTrivial(), rep_r4.getTrivial(), rep_r5.getTrivial());

        const ScalarType rep_r7 = rep_r5 * rep_r2;
        const PacketType rep1 = rep * rep_r7;
        const PacketType const1(-42 * c6, -72 * c8, 90 * c9, -110 * c10);
        const ScalarType d2_g = (rep1 * const1).horizontal_add();
        if (r < cutoff) {
            const ScalarType rep_r6 = rep_r5 * rep_r;
            const PacketType rep2 = rep * rep_r6;
            const PacketType const2(12 * c6, 16 * c8, -18 * c9, 20 * c10);
            const ScalarType d_g = (rep2 * const2).horizontal_add();

            const PacketType rep3 = rep * rep_r5;
            const PacketType const3(-c6, -c8, c9, -c10);
            const ScalarType g = (rep3 * const3).horizontal_add();

            result += d2_g;
            result -= d_g * (rep_r3 * (2 * squaredCutoff) - rep_r2 * (2 * cutoff));
            const PacketType rep4(rep_r3.getTrivial(), rep_r4.getTrivial(), rep_r5.getTrivial(), rep_r6.getTrivial());
            const PacketType const4(4 * cutoff, -2 * squaredCutoff, -8 * cutoff * squaredCutoff, 4 * squaredCutoff * squaredCutoff);
            result += g * (rep4 * const4).horizontal_add();
            result *= exp(-square(rep_r * cutoff - 1));
        }
        else
            result += d2_g;
        return result;
    }
}
