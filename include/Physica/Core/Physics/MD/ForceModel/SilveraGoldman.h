/*
 * Copyright 2023-2024 Weibo He.
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
    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell = false>
    class SilveraGoldman final : public PairModel<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>> {
        constexpr static double alpha = 1.713;
        constexpr static double beta = 1.5671;
        constexpr static double gamma = 0.00993;
        constexpr static double cutoff = 8.32;
        constexpr static double c6 = 12.14;
        constexpr static double c8 = 215.2;
        constexpr static double c9 = 143.1;
        constexpr static double c10 = 4813.9;
        constexpr static double squaredCutoff = cutoff * cutoff;

        using This = SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::ValueType;
    public:
        SilveraGoldman(ValueType cutoff_);
        SilveraGoldman(const SilveraGoldman&) = default;
        SilveraGoldman(SilveraGoldman&&) noexcept = default;
        ~SilveraGoldman() = default;
        /* Operators */
        SilveraGoldman& operator=(SilveraGoldman obj) noexcept;
        /* Operations */
        [[nodiscard]] inline T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] inline T force_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] inline T forceConst_functor(T r, T r2) const;
        inline void swap(SilveraGoldman& __restrict obj) noexcept;
    };

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::SilveraGoldman(ValueType cutoff_) : Base(std::move(cutoff_)) {}

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>&
    SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::operator=(SilveraGoldman obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    inline void SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::swap(SilveraGoldman& __restrict obj) noexcept {
        Base::swap(obj);
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    inline T SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, T r, T r2) const {
        T result = exp(-r2 * ValueType(gamma) - r * ValueType(beta) + ValueType(alpha));
        const T rep_r = reciprocal(r);
        const T rep_r2 = square(rep_r);
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r8 = square(rep_r4);
        const T rep_r9 = rep_r8 * rep_r;
        const T rep_r10 = rep_r6 * rep_r4;
        const T g = rep_r6 * ValueType(c6) + rep_r8 * ValueType(c8) - rep_r9 * ValueType(c9) + rep_r10 * ValueType(c10);

        if (r < ValueType(cutoff)) {
            const T f_cutoff = exp(-square(rep_r * ValueType(cutoff) - ValueType(1)));
            result -= g * f_cutoff;
        }
        else
            result -= g;
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    inline T SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, T r, T r2) const {
        const T factor = r * ValueType(gamma * 2) + ValueType(beta);
        T result = exp(-r2 * ValueType(gamma) - (r * ValueType(beta) - ValueType(alpha))) * factor;
        const T rep_r = reciprocal(r);
        const T rep_r2 = square(rep_r);
        const T rep_r3 = rep_r * rep_r2;
        const T rep_r4 = square(rep_r2);
        const T rep_r5 = rep_r2 * rep_r3;
        const T rep_r6 = rep_r5 * rep_r;

        const Vector4D<T> rep{rep_r, rep_r3, rep_r4, rep_r5};
        const Vector4D<T> rep1 = rep * rep_r6;
        const Vector4D<T> const1{-6 * c6, -8 * c8, 9 * c9, -10 * c10};
        const T d_g = (rep1 * const1).sum();

        if (r < ValueType(cutoff)) {
            const T g = (rep_r2 * ValueType(c6) + rep_r4 * ValueType(c8)) * rep_r4 - (rep_r5 * ValueType(c9) - rep_r6 * ValueType(c10)) * rep_r4;
            const T f_cutoff = exp(-square(rep_r * ValueType(cutoff) - ValueType(1)));
            result += (d_g + g * (rep_r3 * ValueType(2 * squaredCutoff) - rep_r2 * ValueType(2 * cutoff))) * f_cutoff;
        }
        else
            result += d_g;
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    inline T SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::forceConst_functor(T r, T r2) const {
        const T factor = square(r * ValueType(2 * gamma) + ValueType(beta)) - ValueType(2 * gamma);
        T result = exp(-r2 * ValueType(gamma) - (r * ValueType(beta) - ValueType(alpha))) * factor;
        const T rep_r = reciprocal(r);
        const T rep_r2 = square(rep_r);
        const T rep_r3 = rep_r * rep_r2;
        const T rep_r4 = square(rep_r2);
        const T rep_r5 = rep_r2 * rep_r3;
        const Vector4D<T> rep{rep_r, rep_r3, rep_r4, rep_r5};

        const T rep_r7 = rep_r5 * rep_r2;
        const Vector4D<T> rep1 = rep * rep_r7;
        const Vector4D<T> const1{-42 * c6, -72 * c8, 90 * c9, -110 * c10};
        const T d2_g = rep1 * const1;
        if (r < ValueType(cutoff)) {
            const T rep_r6 = rep_r5 * rep_r;
            const Vector4D<T> rep2 = rep * rep_r6;
            const Vector4D<T> const2{-12 * c6, -16 * c8, 18 * c9, -20 * c10};
            const T d_g = rep2 * const2;

            const Vector4D<T> rep3 = rep * rep_r5;
            const Vector4D<T> const3{-c6, -c8, c9, -c10};
            const T g = rep3 * const3;

            const T term2 = -d_g * (rep_r3 * ValueType(2 * squaredCutoff) - rep_r2 * ValueType(2 * cutoff));
            const Vector4D<T> rep4{rep_r3, rep_r4, rep_r5, rep_r6};
            const Vector4D<T> const4{4 * cutoff, -2 * squaredCutoff, -8 * cutoff * squaredCutoff, 4 * squaredCutoff * squaredCutoff};
            const T term3 = g * (rep4 * const4);
            result += (d2_g + term2 + term3) * exp(-square(rep_r * ValueType(cutoff) - ValueType(1)));
        }
        else
            result += d2_g;
        return result;
    }
}

namespace Physica {
    template<Scalar T, bool B1, bool B2>
    class Traits<Core::SilveraGoldman<T, B1, B2>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPeriodBoundary = B1;
        constexpr static bool IsSmallCell = B2;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsContractable = false;
    };
}
