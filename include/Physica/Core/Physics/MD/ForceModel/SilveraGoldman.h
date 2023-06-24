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
    /**
     * Potential that suits para-hydrogen
     * 
     * Reference:
     * [1] I. F. Silvera and V. V. Goldman, J. Chem. Phys. 69, 4209 (1978).
     */
    template<class ScalarType, class PosScalarType>
    class SilveraGoldman final : public PairModel<ScalarType, PosScalarType, SilveraGoldman<ScalarType, PosScalarType>> {
        constexpr static double alpha = 1.713;
        constexpr static double beta = 1.5671;
        constexpr static double gamma = 0.00993;
        constexpr static double cutoff = 8.32;
        constexpr static double c6 = 12.14;
        constexpr static double c8 = 215.2;
        constexpr static double c9 = 143.1;
        constexpr static double c10 = 4813.9;

        using This = SilveraGoldman<ScalarType, PosScalarType>;
        using Base = PairModel<ScalarType, PosScalarType, This>;
    public:
        SilveraGoldman(ScalarType cutoff_);
        SilveraGoldman(const SilveraGoldman&) = default;
        SilveraGoldman(SilveraGoldman&&) noexcept = default;
        ~SilveraGoldman() = default;
        /* Operators */
        SilveraGoldman& operator=(SilveraGoldman obj) noexcept;
        /* Operations */
        void swap(SilveraGoldman& obj) noexcept;
        /* Static members */
        static inline ScalarType force_functor(ScalarType r, ScalarType r2);
        static inline ScalarType pot_functor(ScalarType r, ScalarType r2);
    };

    template<class ScalarType, class PosScalarType>
    SilveraGoldman<ScalarType, PosScalarType>::SilveraGoldman(ScalarType cutoff_) : Base(cutoff_) {}

    template<class ScalarType, class PosScalarType>
    SilveraGoldman<ScalarType, PosScalarType>& SilveraGoldman<ScalarType, PosScalarType>::operator=(SilveraGoldman obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    void SilveraGoldman<ScalarType, PosScalarType>::swap(SilveraGoldman& obj) noexcept {
        Base::swap(obj);
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType SilveraGoldman<ScalarType, PosScalarType>::force_functor(ScalarType r, ScalarType r2) {
        using PacketType = typename Internal::BestPacket<ScalarType, 4>::Type;
        static_assert(!std::is_same<ScalarType, PacketType>::value, "[Error]: SIMD is inavailable, implementation must be revised");

        const ScalarType factor = r * (gamma * 2) + beta;
        ScalarType result = exp(-r2 * gamma - (r * beta - alpha)) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;
        const PacketType rep(rep_r.getTrivial(), rep_r3.getTrivial(), rep_r4.getTrivial(), rep_r5.getTrivial());

        const ScalarType rep_r6 = rep_r5 * rep_r;
        const PacketType rep1 = rep * rep_r6.getTrivial();
        const PacketType const1(-6 * c6, -8 * c8, 9 * c9, -10 * c10);
        const ScalarType d_g = horizontal_add(rep1 * const1);
        if (r < cutoff) {
            const ScalarType g = (rep_r2 * c6 + rep_r4 * c8) * rep_r4 - (rep_r5 * c9 + rep_r6 * c10) * rep_r4;
            const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
            result += (d_g + g * (rep_r3 * (2 * cutoff * cutoff) - rep_r2 * (2 * cutoff))) * f_cutoff;
        }
        else {
            result += d_g;
        }
        return result;
    }
    
    template<class ScalarType, class PosScalarType>
    inline ScalarType SilveraGoldman<ScalarType, PosScalarType>::pot_functor(ScalarType r, ScalarType r2) {
        ScalarType result = exp(-r2 * gamma - r * beta + alpha);
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r8 = square(rep_r4);
        const ScalarType rep_r9 = rep_r8 * rep_r;
        const ScalarType rep_r10 = rep_r6 * rep_r4;
        const ScalarType g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;

        if (r < cutoff) {
            const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
            result -= g * f_cutoff;
        }
        else
            result -= g;
        return result;
    }
}
