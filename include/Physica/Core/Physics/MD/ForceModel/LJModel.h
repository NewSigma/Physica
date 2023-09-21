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
    template<class ScalarType, class PosScalarType> class LJModel;

    namespace Internal {
        template<class T, class U>
        class Traits<LJModel<T, U>> {
        public:
            using ScalarType = T;
            using PosScalarType = U;
            constexpr static bool IsPotDependOnAtomIndex = false;
        };
    }

    template<class ScalarType, class PosScalarType>
    class LJModel : public PairModel<LJModel<ScalarType, PosScalarType>> {
        using This = LJModel<ScalarType, PosScalarType>;
        using Base = PairModel<This>;

        ScalarType sigma;
    public:
        LJModel(ScalarType sigma_, ScalarType cutoff_);
        LJModel(const LJModel&) = default;
        LJModel(LJModel&&) noexcept = default;
        ~LJModel() = default;
        /* Operators */
        LJModel& operator=(LJModel obj) noexcept;
        /* Operations */
        void swap(LJModel& obj) noexcept;
        /* Static members */
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, PosScalarType r, PosScalarType r2) const;
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, PosScalarType r, PosScalarType r2) const;
    };

    template<class ScalarType, class PosScalarType>
    LJModel<ScalarType, PosScalarType>::LJModel(ScalarType sigma_, ScalarType cutoff_) : Base(cutoff_), sigma(sigma_) {}

    template<class ScalarType, class PosScalarType>
    LJModel<ScalarType, PosScalarType>& LJModel<ScalarType, PosScalarType>::operator=(LJModel obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType LJModel<ScalarType, PosScalarType>::force_functor(\
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, PosScalarType r, [[maybe_unused]] PosScalarType r2) const {
        const ScalarType rep_r = ScalarType(sigma) / r;
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r7 = rep_r6 * rep_r;
        const ScalarType rep_r13 = rep_r7 * rep_r6;
        return rep_r13 * 2 - rep_r7;
    }

    template<class ScalarType, class PosScalarType>
    inline ScalarType LJModel<ScalarType, PosScalarType>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, [[maybe_unused]] PosScalarType r, PosScalarType r2) const {
        const ScalarType rep_r2 = ScalarType(sigma * sigma) / r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r12 = square(rep_r6);
        return rep_r12 - rep_r6;
    }

    template<class ScalarType, class PosScalarType>
    void LJModel<ScalarType, PosScalarType>::swap(LJModel& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sigma.swap(obj.sigma);
    }
}
