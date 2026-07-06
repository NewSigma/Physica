/*
 * Copyright 2023-2026 Weibo He.
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

namespace Physica {
    /**
     * Potential that suits para-hydrogen
     *
     * Reference:
     * [1] J. Chem. Phys. 69, 4209 (1978); https://doi.org/10.1063/1.437103
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
        using DVector = CoDiff<Vector4D<T>>;
        using typename Base::Tv;
    public:
        SilveraGoldman(Tv cutoff_);
        SilveraGoldman(const This&) = default;
        SilveraGoldman(This&&) noexcept = default;
        ~SilveraGoldman() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] CoDiff<T> pot_functor(size_t i, size_t j, const T& r, const T& r2) const;
        [[nodiscard]] CoDiff<T> force_functor(size_t i, size_t j, const T& r, const T& r2) const;
        [[nodiscard]] CoDiff<T> forceConst_functor(const T& r, const T& r2) const;
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return IsPeriodBoundary; }
    };

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::SilveraGoldman(Tv cutoff_) : Base(std::move(cutoff_)) {}

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    void SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::swap(This& __restrict obj) noexcept {
        Base::swap(obj);
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    CoDiff<T> SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::pot_functor(size_t, size_t, const T& r, const T& r2) const {
        const auto result0 = exp(-r2 * Tv(gamma) - r * Tv(beta) + Tv(alpha));
        const auto rep_r = reciprocal(r);
        const auto rep_r2 = square(rep_r);
        const auto rep_r4 = square(rep_r2);
        const auto rep_r6 = rep_r4 * rep_r2;
        const auto rep_r8 = square(rep_r4);
        const auto rep_r9 = rep_r8 * rep_r;
        const auto rep_r10 = rep_r6 * rep_r4;
        const auto g = rep_r6 * Tv(c6) + rep_r8 * Tv(c8) - rep_r9 * Tv(c9) + rep_r10 * Tv(c10);

        CoDiff<T> result;
        if (r < Tv(cutoff)) {
            const auto f_cutoff = exp(-square(rep_r * Tv(cutoff) - Tv(1)));
            result = result0 - g * f_cutoff;
        }
        else
            result = result0 - g;
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    CoDiff<T> SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::force_functor(size_t, size_t, const T& r, const T& r2) const {
        const auto factor = r * Tv(gamma * 2) + Tv(beta);
        const auto result0 = exp(-r2 * Tv(gamma) - (r * Tv(beta) - Tv(alpha))) * factor;
        const auto rep_r = reciprocal(r);
        const auto rep_r2 = square(rep_r);
        const auto rep_r3 = rep_r * rep_r2;
        const auto rep_r4 = square(rep_r2);
        const auto rep_r5 = rep_r2 * rep_r3;
        const auto rep_r6 = rep_r5 * rep_r;

        const Vector4D<T> rep{rep_r.value(), rep_r3.value(), rep_r4.value(), rep_r5.value()};
        const DVector rep1 = rep * rep_r6;
        const Vector4D<T> const1{-6 * c6, -8 * c8, 9 * c9, -10 * c10};
        const auto d_g = rep1 * const1;

        CoDiff<T> result;
        if (r < Tv(cutoff)) {
            const auto g = (rep_r2 * Tv(c6) + rep_r4 * Tv(c8)) * rep_r4 - (rep_r5 * Tv(c9) - rep_r6 * Tv(c10)) * rep_r4;
            const auto f_cutoff = exp(-square(rep_r * Tv(cutoff) - Tv(1)));
            result = result0 + (d_g + g * (rep_r3 * Tv(2 * squaredCutoff) - rep_r2 * Tv(2 * cutoff))) * f_cutoff;
        }
        else
            result = result0 + d_g;

        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    CoDiff<T> SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>::forceConst_functor(const T& r, const T& r2) const {
        const auto factor = square(r * Tv(2 * gamma) + Tv(beta)) - Tv(2 * gamma);
        const auto result0 = exp(-r2 * Tv(gamma) - (r * Tv(beta) - Tv(alpha))) * factor;
        const auto rep_r = reciprocal(r);
        const auto rep_r2 = square(rep_r);
        const auto rep_r3 = rep_r * rep_r2;
        const auto rep_r4 = square(rep_r2);
        const auto rep_r5 = rep_r2 * rep_r3;
        const Vector4D<T> rep{rep_r.value(), rep_r3.value(), rep_r4.value(), rep_r5.value()};

        const auto rep_r7 = rep_r5 * rep_r2;
        const DVector rep1 = rep * rep_r7;
        const Vector4D<Tv> const1{-42 * c6, -72 * c8, 90 * c9, -110 * c10};
        const auto d2_g = rep1 * const1;

        CoDiff<T> result;
        if (r < Tv(cutoff)) {
            const auto rep_r6 = rep_r5 * rep_r;
            const DVector rep2 = rep * rep_r6;
            const Vector4D<Tv> const2{-12 * c6, -16 * c8, 18 * c9, -20 * c10};
            const auto d_g = rep2 * const2;

            const DVector rep3 = rep * rep_r5;
            const Vector4D<Tv> const3{-c6, -c8, c9, -c10};
            const auto g = rep3 * const3;

            const auto term2 = -d_g * (rep_r3 * Tv(2 * squaredCutoff) - rep_r2 * Tv(2 * cutoff));
            const Vector4D<T> rep4{rep_r3.value(), rep_r4.value(), rep_r5.value(), rep_r6.value()};
            const Vector4D<Tv> const4{4 * cutoff, -2 * squaredCutoff, -8 * cutoff * squaredCutoff, 4 * squaredCutoff * squaredCutoff};
            const auto term3 = g * (rep4 * const4);
            result = result0 + (d2_g + term2 + term3) * exp(-square(rep_r * Tv(cutoff) - Tv(1)));
        }
        else
            result = result0 + d2_g;

        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }
}

namespace Physica {
    template<Scalar T, bool B1, bool B2>
    class Traits<SilveraGoldman<T, B1, B2>> {
    public:
        using ScalarType = T;
        constexpr static bool IsSmallCell = B2;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsContractable = false;
    };
}
