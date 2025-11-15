/*
 * Copyright 2025 Weibo He.
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

#include "DQMCImpl/ImagKinetic.h"
#include "DQMCImpl/GreenProd.h"

namespace Physica {
    /**
     * Reference:
     * [1] Phys. Rev. B 40, 506; https://doi.org/10.1103/PhysRevB.40.506
     */
    template<Scalar T>
    class DQMC {
        using This = DQMC<T>;
        using Params = HubbardParams<T>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        constexpr static bool isComplex = T::isComplex;
    private:
        const Params* params;
        ImagKinetic<T> kinetic;
        GreenProd<T> productor;

        VectorND<Tr> probs;
        Array<int> sites;
        int cursor = 0;
    public:
        DQMC() = delete;
        DQMC(const Params& params_);
        DQMC(const This&) = default;
        DQMC(This&&) noexcept = default;
        ~DQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R = Random<>>
        void step_random();

        template<RNG R = Random<>>
        void step_spin();
        template<RNG R = Random<>>
        void step_spin_for(int numStep);

        void calcGreens(int split);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return *params; }
        [[nodiscard]] int getNumSite() const noexcept { return kinetic.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return kinetic.getNumSplit(); }
        [[nodiscard]] const auto& getAuxField() const noexcept { return kinetic.getAuxField(); }
        [[nodiscard]] auto& getGreens() noexcept { return kinetic.getGreens(); }
        [[nodiscard]] Tv getSign() const noexcept { return kinetic.getSign(); }
    private:
        /* Operations */
        void metropolis(int site, int split, Tr prob);
        /* Static members */
        template<RNG R>
        [[nodiscard]] static Array<int> makeRandomSites(int numSite);
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_)
            : params(&params_)
            , kinetic(params_.getNumSite(), params_.getNumSplit())
            , productor(params_)
            , probs(params_.getNumSite())
            , sites(params_.getNumSite()) {
        assert(getNumSplit() % 2 == 0 && "[Error]: An even number of splits is required");
        for (int i = 0; i < getNumSite(); ++i)
            sites[i] = i;
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_random() {
        kinetic.template random_uniform<R>();
        productor.invalidates(getAuxField(), params->getAlpha());
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin() {
        std::ranges::shuffle(sites, R::getInstance());
        probs.template random_uniform<R>();
        calcGreens(cursor);
        for (int i = 0; i < getNumSite(); ++i)
            metropolis(sites[i], cursor, probs[i]);
        productor.invalidate(cursor);
        cursor = (cursor + 1) % getNumSplit();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int _ = 0; _ < numStep; ++_)
            step_spin<R>();
    }

    template<Scalar T>
    void DQMC<T>::calcGreens(int split) {
        productor.calcGreens(getGreens(), split, params->calcBetaMu());
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        kinetic.swap(obj.kinetic);
        productor.swap(obj.productor);
        sites.swap(obj.sites);
        std::swap(cursor, obj.cursor);
    }

    template<Scalar T>
    void DQMC<T>::metropolis(int site, int split, Tr prob) {
        const bool accept = prob < kinetic.calcP(site, split, params->getAlpha());
        if (accept) {
            const Tr x = Tr(2) * params->getAlpha() * getAuxField()(site, split);
            const Vector2D<Tr> arr = exp(Vector2D<Tr>{-x, x});

            auto deltas = kinetic.calcDelta(site, split, params->getAlpha());
            auto ratios = kinetic.calcRatio(site, deltas);
            productor.single_flip(site, split, arr, ratios);
            kinetic.single_flip(site, split, params->getAlpha());
        }
    }
}
