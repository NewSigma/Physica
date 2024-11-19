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
    template<Scalar T, bool IsSmallCell = false>
    class LJModel : public PairModel<LJModel<T, IsSmallCell>> {
        using This = LJModel<T, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::ValueType;

        T sigma;
        T sigma1;
    public:
        LJModel(T sigma_, ValueType cutoff);
        LJModel(const LJModel&) = default;
        LJModel(LJModel&&) noexcept = default;
        ~LJModel() = default;
        /* Operators */
        LJModel& operator=(LJModel obj) noexcept;
        /* Operations */
        void swap(LJModel& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] inline T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] inline T force_functor(size_t i, size_t j, T r, T r2) const;
    };

    template<Scalar T, bool IsSmallCell>
    LJModel<T, IsSmallCell>::LJModel(T sigma_, ValueType cutoff)
            : Base(), sigma(std::move(sigma_)) {
        sigma1 = ValueType(6) / sigma;
        Base::setCutoff(std::move(cutoff));
    }

    template<Scalar T, bool IsSmallCell>
    LJModel<T, IsSmallCell>& LJModel<T, IsSmallCell>::operator=(LJModel obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, bool IsSmallCell>
    inline T LJModel<T, IsSmallCell>::pot_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, [[maybe_unused]] T r, T r2) const {
        const T rep_r2 = T(sigma * sigma) / r2;
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r12 = square(rep_r6);
        return rep_r12 - rep_r6;
    }

    template<Scalar T, bool IsSmallCell>
    inline T LJModel<T, IsSmallCell>::force_functor(
            [[maybe_unused]] size_t i, [[maybe_unused]] size_t j, T r, [[maybe_unused]] T r2) const {
        const T rep_r = sigma / r;
        const T rep_r2 = square(rep_r);
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r7 = rep_r6 * rep_r;
        const T rep_r13 = rep_r7 * rep_r6;
        return (rep_r13 * ValueType(2) - rep_r7) * sigma1;
    }

    template<Scalar T, bool IsSmallCell>
    void LJModel<T, IsSmallCell>::swap(LJModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sigma.swap(obj.sigma);
        sigma1.swap(obj.sigma1);
    }
}

namespace Physica {
    template<class T, bool B>
    class Traits<LJModel<T, B>> : public Traits<PairModel<LJModel<T, B>>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsSmallCell = B;
        constexpr static bool IsContractable = false;
    };
}
