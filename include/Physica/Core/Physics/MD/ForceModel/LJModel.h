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
    template<Scalar T, bool IsSmallCell = false>
    class LJModel : public PairModel<LJModel<T, IsSmallCell>> {
        using This = LJModel<T, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::Tv;

        T sigma;
        CoDiff<T> sigma1;
    public:
        LJModel(T sigma_, Tv cutoff);
        LJModel(const This&) = default;
        LJModel(This&&) noexcept = default;
        ~LJModel() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] CoDiff<T> pot_functor(size_t i, size_t j, const T& r, const T& r2) const;
        [[nodiscard]] CoDiff<T> force_functor(size_t i, size_t j, const T& r, const T& r2) const;
    };

    template<Scalar T, bool IsSmallCell>
    LJModel<T, IsSmallCell>::LJModel(T sigma_, Tv cutoff)
            : Base(), sigma(std::move(sigma_)) {
        sigma1 = Tv(6) / sigma;
        Base::setCutoff(std::move(cutoff));
    }

    template<Scalar T, bool IsSmallCell>
    CoDiff<T> LJModel<T, IsSmallCell>::pot_functor(size_t, size_t, const T&, const T& r2) const {
        auto rep_r2 = square(sigma) / r2;
        auto rep_r4 = square(rep_r2);
        auto rep_r6 = rep_r4 * rep_r2;
        auto rep_r12 = square(rep_r6);
        auto result = rep_r12 - rep_r6;
        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }

    template<Scalar T, bool IsSmallCell>
    CoDiff<T> LJModel<T, IsSmallCell>::force_functor(size_t, size_t, const T& r, const T&) const {
        auto rep_r = sigma / r;
        auto rep_r2 = square(rep_r);
        auto rep_r4 = square(rep_r2);
        auto rep_r6 = rep_r4 * rep_r2;
        auto rep_r7 = rep_r6 * rep_r;
        auto rep_r13 = rep_r7 * rep_r6;
        auto result = (rep_r13 * Tv(2) - rep_r7) * sigma1;
        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }

    template<Scalar T, bool IsSmallCell>
    void LJModel<T, IsSmallCell>::swap(This& __restrict obj) noexcept {
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
