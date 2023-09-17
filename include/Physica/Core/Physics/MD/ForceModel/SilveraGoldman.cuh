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

#include "PairModel.cuh"
#include "SilveraGoldman.h"

namespace Physica::Core {
    namespace Internal {
        template<class T, class U>
        class Traits<device_obj<SilveraGoldman<T, U>>> {
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
    class device_obj<SilveraGoldman<ScalarType, PosScalarType>> final : public device_obj<PairModel<SilveraGoldman<ScalarType, PosScalarType>>> {
        constexpr static float alpha = 1.713;
        constexpr static float beta = 1.5671;
        constexpr static float gamma = 0.00993;
        constexpr static float cutoff = 8.32;
        constexpr static float c6 = 12.14;
        constexpr static float c8 = 215.2;
        constexpr static float c9 = 143.1;
        constexpr static float c10 = 4813.9;

        using host_obj = SilveraGoldman<ScalarType, PosScalarType>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<PairModel<host_obj>>;
    public:
        device_obj(size_t numParticle, ScalarType cutoff_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] __host__ __device__ inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] __host__ __device__ inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        void swap(device_obj& obj) noexcept;
    };

    template<class ScalarType, class PosScalarType>
    device_obj<SilveraGoldman<ScalarType, PosScalarType>>::device_obj(size_t numParticle, ScalarType cutoff_)
            : Base(numParticle, cutoff_) {}

    template<class ScalarType, class PosScalarType>
    void device_obj<SilveraGoldman<ScalarType, PosScalarType>>::swap(device_obj& obj) noexcept {
        Base::swap(obj);
    }

    template<class ScalarType, class PosScalarType>
    __host__ __device__ inline ScalarType device_obj<SilveraGoldman<ScalarType, PosScalarType>>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, ScalarType r2) const {
        const ScalarType factor = r * (gamma * 2) + beta;
        ScalarType result = exp(-r2 * gamma - (r * beta - alpha)) * factor;
        const ScalarType rep_r = reciprocal(r);
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r3 = rep_r * rep_r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r5 = rep_r2 * rep_r3;

        const ScalarType rep_r6 = rep_r5 * rep_r;
        ScalarType d_g = ScalarType(-6 * c6) * rep_r + ScalarType(-8 * c8) * rep_r3 + ScalarType(9 * c9) * rep_r4 + ScalarType(-10 * c10) * rep_r5;
        d_g *= rep_r6;
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
    __host__ __device__ inline ScalarType device_obj<SilveraGoldman<ScalarType, PosScalarType>>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, ScalarType r2) const {
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
