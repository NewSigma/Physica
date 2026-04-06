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

#include "PairModel.cuh"
#include "SilveraGoldman.h"

namespace Physica {
    /**
     * Potential that suits para-hydrogen
     * 
     * Reference:
     * [1] J. Chem. Phys. 69, 4209 (1978); https://doi.org/10.1063/1.437103
     */
    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    class device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>> final
            : public device_obj<PairModel<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>> {
        constexpr static float alpha = 1.713;
        constexpr static float beta = 1.5671;
        constexpr static float gamma = 0.00993;
        constexpr static float cutoff = 8.32;
        constexpr static float c6 = 12.14;
        constexpr static float c8 = 215.2;
        constexpr static float c9 = 143.1;
        constexpr static float c10 = 4813.9;

        using host_obj = SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<PairModel<host_obj>>;
        using typename Base::MDCellType;
    public:
        device_obj() = default;
        device_obj(size_t numParticle, T cutoff_);
        device_obj(const host_obj& obj, size_t numParticle);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] __host__ __device__ T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] __host__ __device__ T force_functor(size_t i, size_t j, T r, T r2) const;
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return IsPeriodBoundary; }
    };

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>::device_obj(size_t numParticle, T cutoff_)
            : Base(numParticle, cutoff_) {}

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>::device_obj(const host_obj& obj, size_t numParticle)
            : device_obj(numParticle, obj.getCutoff()) {}

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    void device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>::swap(This& __restrict obj) noexcept {
        Base::swap(obj);
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    __host__ __device__ T device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, T r, T r2) const {
        T result = exp(-r2 * gamma - r * beta + alpha);
        const T rep_r = reciprocal(r);
        const T rep_r2 = square(rep_r);
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r8 = square(rep_r4);
        const T rep_r9 = rep_r8 * rep_r;
        const T rep_r10 = rep_r6 * rep_r4;
        const T g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;

        if (r < cutoff) {
            const T f_cutoff = exp(-square(rep_r * cutoff - 1));
            result -= g * f_cutoff;
        }
        else
            result -= g;
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary, bool IsSmallCell>
    __host__ __device__ T device_obj<SilveraGoldman<T, IsPeriodBoundary, IsSmallCell>>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, T r, T r2) const {
        const T factor = r * (gamma * 2) + beta;
        T result = exp(-r2 * gamma - (r * beta - alpha)) * factor;
        const T rep_r = reciprocal(r);
        const T rep_r2 = square(rep_r);
        const T rep_r3 = rep_r * rep_r2;
        const T rep_r4 = square(rep_r2);
        const T rep_r5 = rep_r2 * rep_r3;

        const T rep_r6 = rep_r5 * rep_r;
        T d_g = T(-6 * c6) * rep_r + T(-8 * c8) * rep_r3 + T(9 * c9) * rep_r4 + T(-10 * c10) * rep_r5;
        d_g *= rep_r6;
        if (r < cutoff) {
            const T g = (rep_r2 * c6 + rep_r4 * c8) * rep_r4 - (rep_r5 * c9 - rep_r6 * c10) * rep_r4;
            const T f_cutoff = exp(-square(rep_r * cutoff - 1));
            result += (d_g + g * (rep_r3 * (2 * cutoff * cutoff) - rep_r2 * (2 * cutoff))) * f_cutoff;
        }
        else {
            result += d_g;
        }
        return result;
    }
}

namespace Physica {
    template<class T, bool B, bool IsSmallCell>
    class Traits<device_obj<SilveraGoldman<T, B, IsSmallCell>>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsContractable = false;
    };
}
