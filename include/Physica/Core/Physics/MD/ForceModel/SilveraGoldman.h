/*
 * Copyright 2023-2024 WeiBo He.
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
    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell = false>
    class SilveraGoldman final : public PairModel<SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>> {
        constexpr static double alpha = 1.713;
        constexpr static double beta = 1.5671;
        constexpr static double gamma = 0.00993;
        constexpr static double cutoff = 8.32;
        constexpr static double c6 = 12.14;
        constexpr static double c8 = 215.2;
        constexpr static double c9 = 143.1;
        constexpr static double c10 = 4813.9;
        constexpr static double squaredCutoff = cutoff * cutoff;

        using This = SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::PlainScalar;
        using Vector4D = Vector<ScalarType, 4>;
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
        inline void swap(SilveraGoldman& __restrict obj) noexcept;
    };

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::SilveraGoldman(PlainScalar cutoff_) : Base(std::move(cutoff_)) {}

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>&
    SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::operator=(SilveraGoldman obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    inline void SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::swap(SilveraGoldman& __restrict obj) noexcept {
        Base::swap(obj);
    }

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    inline ScalarType SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::pot_functor(
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

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    inline ScalarType SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, ScalarType r2) const {
        const ScalarType factor = r * PlainScalar(gamma * 2) + PlainScalar(beta);
        ScalarType result = exp(-r2 * PlainScalar(gamma) - (r * PlainScalar(beta) - PlainScalar(alpha))) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;
        const ScalarType rep_r6 = rep_r5 * rep_r;

        const Vector4D rep{rep_r, rep_r3, rep_r4, rep_r5};
        const Vector4D rep1 = rep * rep_r6;
        const Vector4D const1{-6 * c6, -8 * c8, 9 * c9, -10 * c10};
        const ScalarType d_g = (rep1 * const1).horizontal_add();

        if (r < PlainScalar(cutoff)) {
            const ScalarType g = (rep_r2 * PlainScalar(c6) + rep_r4 * PlainScalar(c8)) * rep_r4 - (rep_r5 * PlainScalar(c9) - rep_r6 * PlainScalar(c10)) * rep_r4;
            const ScalarType f_cutoff = exp(-square(rep_r * PlainScalar(cutoff) - PlainScalar(1)));
            result += (d_g + g * (rep_r3 * PlainScalar(2 * squaredCutoff) - rep_r2 * PlainScalar(2 * cutoff))) * f_cutoff;
        }
        else
            result += d_g;
        return result;
    }

    template<class ScalarType, bool IsPeriodBoundary, bool IsSmallCell>
    inline ScalarType SilveraGoldman<ScalarType, IsPeriodBoundary, IsSmallCell>::forceConst_functor(ScalarType r, ScalarType r2) const {
        const ScalarType factor = square(r * PlainScalar(2 * gamma) + PlainScalar(beta)) - PlainScalar(2 * gamma);
        ScalarType result = exp(-r2 * PlainScalar(gamma) - (r * PlainScalar(beta) - PlainScalar(alpha))) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;
        const Vector4D rep{rep_r, rep_r3, rep_r4, rep_r5};

        const ScalarType rep_r7 = rep_r5 * rep_r2;
        const Vector4D rep1 = rep * rep_r7;
        const Vector4D const1{-42 * c6, -72 * c8, 90 * c9, -110 * c10};
        const ScalarType d2_g = rep1 * const1;
        if (r < PlainScalar(cutoff)) {
            const ScalarType rep_r6 = rep_r5 * rep_r;
            const Vector4D rep2 = rep * rep_r6;
            const Vector4D const2{-12 * c6, -16 * c8, 18 * c9, -20 * c10};
            const ScalarType d_g = rep2 * const2;

            const Vector4D rep3 = rep * rep_r5;
            const Vector4D const3{-c6, -c8, c9, -c10};
            const ScalarType g = rep3 * const3;

            const ScalarType term2 = -d_g * (rep_r3 * PlainScalar(2 * squaredCutoff) - rep_r2 * PlainScalar(2 * cutoff));
            const Vector4D rep4{rep_r3, rep_r4, rep_r5, rep_r6};
            const Vector4D const4{4 * cutoff, -2 * squaredCutoff, -8 * cutoff * squaredCutoff, 4 * squaredCutoff * squaredCutoff};
            const ScalarType term3 = g * (rep4 * const4);
            result += (d2_g + term2 + term3) * exp(-square(rep_r * PlainScalar(cutoff) - PlainScalar(1)));
        }
        else
            result += d2_g;
        return result;
    }
}

namespace Physica {
    template<class T, bool B1, bool B2>
    class Traits<Core::SilveraGoldman<T, B1, B2>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPeriodBoundary = B1;
        constexpr static bool IsSmallCell = B2;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsContractable = false;
    };
}
