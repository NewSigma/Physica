/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "PairModel.h"

namespace Physica {
    /**
     * \class LJModel1 is a variation of \class LJModel
     *
     * Reference:
     * [1] Phys. Rev. 188, 1407; https://doi.org/10.1103/PhysRev.188.1407
     */
    template<Scalar T, bool IsSmallCell = false>
    class LJModel1 : public PairModel<LJModel1<T, IsSmallCell>> {
        using This = LJModel1<T, IsSmallCell>;
        using Base = PairModel<This>;
        using typename Base::Tv;

        T sigma;
        T epsilon;
        T factor;
    public:
        LJModel1(T sigma_, T epsilon_, Tv cutoff);
        LJModel1(const This&) = default;
        LJModel1(This&&) noexcept = default;
        ~LJModel1() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] T force_functor(size_t i, size_t j, T r, T r2) const;
    };

    template<Scalar T, bool IsSmallCell>
    LJModel1<T, IsSmallCell>::LJModel1(T sigma_, T epsilon_, Tv cutoff)
            : Base(), sigma(std::move(sigma_)), epsilon(std::move(epsilon_)), factor(Tv(12) * epsilon / sigma) {
        Base::setCutoff(std::move(cutoff));
    }

    template<Scalar T, bool IsSmallCell>
    auto LJModel1<T, IsSmallCell>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T, bool IsSmallCell>
    T LJModel1<T, IsSmallCell>::pot_functor(size_t, size_t, T, T r2) const {
        const T rep_r2 = T(sigma * sigma) / r2;
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r12 = square(rep_r6);
        return (rep_r12 - rep_r6 * T(2)) * epsilon;
    }

    template<Scalar T, bool IsSmallCell>
    T LJModel1<T, IsSmallCell>::force_functor(size_t, size_t, T r, T) const {
        const T rep_r = sigma / r;
        const T rep_r2 = square(rep_r);
        const T rep_r4 = square(rep_r2);
        const T rep_r6 = rep_r4 * rep_r2;
        const T rep_r7 = rep_r6 * rep_r;
        const T rep_r13 = rep_r7 * rep_r6;
        return (rep_r13 - rep_r7) * factor;
    }

    template<Scalar T, bool IsSmallCell>
    void LJModel1<T, IsSmallCell>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sigma.swap(obj.sigma);
        epsilon.swap(obj.epsilon);
        factor.swap(obj.factor);
    }
}

namespace Physica {
    template<class T, bool B>
    class Traits<LJModel1<T, B>> : public Traits<PairModel<LJModel1<T, B>>> {
    public:
        using ScalarType = T;
        constexpr static double IsSmallCell = B;
        constexpr static bool IsPotDependOnAtomIndex = false;
        constexpr static bool IsLatticeDependent = false;
        constexpr static bool IsContractable = false;
    };
}
