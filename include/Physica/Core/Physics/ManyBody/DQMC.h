/*
 * Copyright 2025-2026 Weibo He.
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
        using MaybeScore = std::conditional_t<T::isDiffable(), Trv, Empty>;

        constexpr static bool isComplex = T::isComplex();
        static_assert(T::Prec == Float64, "[Warn]: It is highly recommended to use high-precision floats");
    private:
        const Params* params;
        ImagKinetic<T> kinetic;
        GreenProd<T> productor;

        VectorND<Trv> probs;
        Array<int> sites;
        int cursor = 0;

        [[no_unique_address]] Trv score = 0;
        uint64_t numTotal = 0;
        uint64_t numAccept = 0;
    public:
        DQMC() = delete;
        explicit DQMC(const Params& params_);
        DQMC(const This&) = default;
        DQMC(This&&) noexcept = default;
        ~DQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        auto step_random();

        template<RNG R, ExecutePolicy P = Sequential>
        void step();
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep);

        template<ExecutePolicy P = Sequential>
        auto calcGreens(int split);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return *params; }
        [[nodiscard]] int getNumSite() const noexcept { return kinetic.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return kinetic.getNumSplit(); }
        [[nodiscard]] const auto& getAuxField() const noexcept { return kinetic.getAuxField(); }
        [[nodiscard]] const auto& getGreens() noexcept { return kinetic.getGreens(); }
        [[nodiscard]] T getRSign() const noexcept;
        [[nodiscard]] uint64_t getNumTotal() const noexcept { return numTotal; }
        [[nodiscard]] uint64_t getNumAccept() const noexcept { return numAccept; }
    private:
        /* Operations */
        void metropolis(int site, int split, Trv prob);
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
    auto DQMC<T>::step_random() {
        kinetic.template random_uniform<R>();
        productor.invalidates(getAuxField(), params->getAlpha());
        cursor = 0;
        return calcGreens(0);
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void DQMC<T>::step() {
        std::ranges::shuffle(sites, R::getInstance());
        probs.template random_uniform<R>();
        for (int i = 0; i < getNumSite(); ++i)
            metropolis(sites[i], cursor, probs[i]);
        productor.invalidate(cursor);
        cursor = (cursor + 1) % getNumSplit();

        calcGreens<P>(cursor);
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void DQMC<T>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int _ = 0; _ < numStep; ++_)
            step<R, P>();

        numTotal = 0;
        numAccept = 0;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto DQMC<T>::calcGreens(int split) {
        if constexpr (T::isDiffable()) {
            const auto [lnAbsDet, _] = productor.template calcDetGreens<P>(kinetic.getGreens(), split, params->calcBetaMu());
            score = lnAbsDet.grad();
        }
        else
            return productor.template calcGreens<P>(kinetic.getGreens(), split, params->calcBetaMu()).wait();
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        kinetic.swap(obj.kinetic);
        productor.swap(obj.productor);
        sites.swap(obj.sites);
        std::swap(score, obj.score);
        std::swap(cursor, obj.cursor);
    }

    template<Scalar T>
    T DQMC<T>::getRSign() const noexcept {
        if constexpr (T::isDiffable()) {
            const Tv rsign = kinetic.getRSign();
            return T(rsign, rsign * score);
        }
        else
            return kinetic.getRSign();
    }

    template<Scalar T>
    void DQMC<T>::metropolis(int site, int split, Trv prob) {
        const bool accept = prob < kinetic.calcP(site, split, params->getAlpha());
        if (accept) {
            const Tr x = Tr(2) * params->getAlpha() * getAuxField()[site, split];
            const Vector2D<Tr> arr = exp(Vector2D<Tr>{-x, x});
            productor.single_flip(site, split, arr);
            kinetic.single_flip(site, split, params->getAlpha());
            numAccept += 1;
        }
        numTotal += 1;
    }
}
