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
    template<class ScalarType, bool IsSmallCell = false> class LJModel1;

    namespace Internal {
        template<class T, bool B>
        class Traits<LJModel1<T, B>> : public Traits<PairModel<LJModel1<T, B>>> {
        public:
            using ScalarType = T;
            constexpr static double IsSmallCell = B;
            constexpr static bool IsPotDependOnAtomIndex = false;
            constexpr static bool IsLatticeDependent = false;
        };
    }
    /**
     * \class LJModel1 is a variation of \class LJModel
     * 
     * Reference:
     * [1] Phys. Rev. 188, 1407; https://doi.org/10.1103/PhysRev.188.1407
     */
    template<class ScalarType, bool IsSmallCell>
    class LJModel1 : public PairModel<LJModel1<ScalarType, IsSmallCell>> {
        using This = LJModel1<ScalarType, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::PlainScalar;

        ScalarType sigma;
        ScalarType epsilon;
        ScalarType factor;
    public:
        LJModel1(ScalarType sigma_, ScalarType epsilon_, PlainScalar cutoff);
        LJModel1(const LJModel1&) = default;
        LJModel1(LJModel1&&) noexcept = default;
        ~LJModel1() = default;
        /* Operators */
        LJModel1& operator=(LJModel1 obj) noexcept;
        /* Operations */
        void swap(LJModel1& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    };

    template<class ScalarType, bool IsSmallCell>
    LJModel1<ScalarType, IsSmallCell>::LJModel1(ScalarType sigma_, ScalarType epsilon_, PlainScalar cutoff)
            : Base(), sigma(std::move(sigma_)), epsilon(std::move(epsilon_)) {
        factor = PlainScalar(12) * epsilon / sigma;
        Base::setCutoff(std::move(cutoff));
    }

    template<class ScalarType, bool IsSmallCell>
    LJModel1<ScalarType, IsSmallCell>& LJModel1<ScalarType, IsSmallCell>::operator=(LJModel1 obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType LJModel1<ScalarType, IsSmallCell>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, [[maybe_unused]] ScalarType r, ScalarType r2) const {
        const ScalarType rep_r2 = ScalarType(sigma * sigma) / r2;
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r12 = square(rep_r6);
        return (rep_r12 - rep_r6 * ScalarType(2)) * epsilon;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType LJModel1<ScalarType, IsSmallCell>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType rep_r = sigma / r;
        const ScalarType rep_r2 = square(rep_r);
        const ScalarType rep_r4 = square(rep_r2);
        const ScalarType rep_r6 = rep_r4 * rep_r2;
        const ScalarType rep_r7 = rep_r6 * rep_r;
        const ScalarType rep_r13 = rep_r7 * rep_r6;
        return (rep_r13 - rep_r7) * factor;
    }

    template<class ScalarType, bool IsSmallCell>
    void LJModel1<ScalarType, IsSmallCell>::swap(LJModel1& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sigma.swap(obj.sigma);
        epsilon.swap(obj.epsilon);
        factor.swap(obj.factor);
    }
}
