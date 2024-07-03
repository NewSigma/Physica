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
    template<class ScalarType, bool IsSmallCell = false>
    class LJModel : public PairModel<LJModel<ScalarType, IsSmallCell>> {
        using This = LJModel<ScalarType, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::PlainScalar;

        ScalarType sigma;
        ScalarType sigma1;
    public:
        LJModel(ScalarType sigma_, PlainScalar cutoff);
        LJModel(const LJModel&) = default;
        LJModel(LJModel&&) noexcept = default;
        ~LJModel() = default;
        /* Operators */
        LJModel& operator=(LJModel obj) noexcept;
        /* Operations */
        void swap(LJModel& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    };

    template<class ScalarType, bool IsSmallCell>
    LJModel<ScalarType, IsSmallCell>::LJModel(ScalarType sigma_, PlainScalar cutoff)
            : Base(), sigma(std::move(sigma_)) {
        sigma1 = PlainScalar(6) / sigma;
        Base::setCutoff(std::move(cutoff));
    }

    template<class ScalarType, bool IsSmallCell>
    LJModel<ScalarType, IsSmallCell>& LJModel<ScalarType, IsSmallCell>::operator=(LJModel obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType LJModel<ScalarType, IsSmallCell>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, [[maybe_unused]] ScalarType r, ScalarType r2) const {
        const ScalarType rep_r2 = ScalarType(sigma * sigma) / r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r12 = square(rep_r6);
        return rep_r12 - rep_r6;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType LJModel<ScalarType, IsSmallCell>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType rep_r = sigma / r;
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r7 = rep_r6 * rep_r;
        const ScalarType rep_r13 = rep_r7 * rep_r6;
        return (rep_r13 * PlainScalar(2) - rep_r7) * sigma1;
    }

    template<class ScalarType, bool IsSmallCell>
    void LJModel<ScalarType, IsSmallCell>::swap(LJModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sigma.swap(obj.sigma);
        sigma1.swap(obj.sigma1);
    }
}

namespace Physica {
    using namespace Core;

    template<class T, bool B>
    class Traits<LJModel<T, B>> : public Traits<PairModel<LJModel<T, B>>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsSmallCell = B;
        constexpr static bool IsContractable = false;
    };
}
