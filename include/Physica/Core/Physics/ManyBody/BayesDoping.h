/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "GreenSampler/ScalarSampler.h"
#include "Physica/Core/ML/BayesOpt.h"
#include "Physica/Core/ML/Optimizer/Adadelta.h"

namespace Physica {
    /**
     * \class BayesDoping: Bayesian optimization for doping
     */
    template<Scalar T>
    class BayesDoping {
        using This = BayesDoping<T>;
        using Tv = T::ValueType;

        using GlobalOpt = Vegas<Tv, true>;
    private:
        GlobalOpt vegas;
        int numBayesIter;
        int numVegasIter;
        int numWarmupDQMC;
        int numSampleDQMC;
        int numStepSGD;

        GPModel<Diff<T, DiffMode::Reverse>> model;
        Adadelta<T> opt;
        VectorND<T> mus;
        VectorND<T> rhos;
        VectorND<T> noises;
    public:
        BayesDoping(GlobalOpt vegas_, int numBayesIter, int numVegasIter, int numWarmupDQMC, int numSampleDQMC, int numStepSGD);
        BayesDoping(const This&) = default;
        BayesDoping(This&&) noexcept = default;
        ~BayesDoping() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<RNG R>
        [[nodiscard]] T solve(T target, HubbardParams<T>& params, auto& dqmc, auto&&... args);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getModel() const noexcept { return model; }
        [[nodiscard]] const auto& getChemMus() const noexcept { return mus; }
        [[nodiscard]] const auto& getDensity() const noexcept { return rhos; }
        [[nodiscard]] const auto& getNoises() const noexcept { return noises; }
        [[nodiscard]] size_t getNumSamples() const noexcept { return mus.getLength(); }
    private:
        static T likelyhood(T x, T y) noexcept { return -square(x - y); }
    };

    template<Scalar T>
    BayesDoping<T>::BayesDoping(GlobalOpt vegas_, int numBayesIter, int numVegasIter, int numWarmupDQMC, int numSampleDQMC, int numStepSGD)
            : vegas(std::move(vegas_))
            , numBayesIter(numBayesIter)
            , numVegasIter(numVegasIter)
            , numWarmupDQMC(numWarmupDQMC)
            , numSampleDQMC(numSampleDQMC)
            , numStepSGD(numStepSGD)
            , model(1) {
        assert(vegas.getDim() == 1 && "[Error]: Doping problem is 1-dim");
        assert(numSampleDQMC > 0 && "[Error]: Must have sample");
        mus.reserve(numBayesIter + 1);
        rhos.reserve(numBayesIter + 1);
        noises.reserve(numBayesIter + 1);

        mus.append(0);
        rhos.append(T(1));
        noises.append(0);
    }

    template<Scalar T>
    template<RNG R>
    auto BayesDoping<T>::solve(T target, HubbardParams<T>& params, auto& dqmc, auto&&... args) -> T {
        assert(T(0) <= target && target <= T(2) && "[Error]: Invalid density");
        const auto runDQMC = [&](T mu) {
            params.setChemMu(mu + params.getRepelU() * 0.5);
            dqmc.template step_random<R>();
            for (int i = 0; i < numWarmupDQMC; ++i)
                dqmc.template step<R>(std::forward<decltype(args)>(args)...);

            ScalarSampler<T> sampler(params, numSampleDQMC);
            for (int i = 0; i < numSampleDQMC; ++i) {
                dqmc.template step<R>(std::forward<decltype(args)>(args)...);
                sampler.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::Density);
            }
            const T rho = sampler.calcMean();
            const T devia = sampler.getObserves().deviation() / sqrt(T(numSampleDQMC));
            return std::pair<T, T>{rho, devia};
        };

        auto likelyhoods = VectorND<T>::generate([this, target](size_t i) {
            return likelyhood(rhos[i], target);
        }, mus.getLength());

        const size_t argmax = likelyhoods.argmax();
        BayesOpt<T> bayes(GlobalOpt(vegas), {mus[argmax]}, likelyhoods[argmax]);
        for (int i = 0; i < numBayesIter; ++i) {
            for (int _ = 0; _ < numStepSGD; ++_) {
                model.regression(mus.transpose(), likelyhoods, square(noises)).reverse(-1);
                model.step(opt);
                model.zero_grad();
            }

            auto [newMu, newRho] = bayes.template propose<R>([&](const VectorND<T>& x) -> T {
                const auto [rho, devia] = runDQMC(x.front());
                auto l = likelyhood(rho, target);
                mus.append(x.front());
                rhos.append(rho);
                noises.append(devia);
                likelyhoods.append(l);
                return l;
            }, model, numVegasIter);
        }
        return bayes.getOptimal().front();
    }

    template<Scalar T>
    void BayesDoping<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        vegas.swap(obj.vegas);
        std::swap(numBayesIter, obj.numBayesIter);
        std::swap(numVegasIter, obj.numVegasIter);
        std::swap(numWarmupDQMC, obj.numWarmupDQMC);
        std::swap(numSampleDQMC, obj.numSampleDQMC);
        std::swap(numStepSGD, obj.numStepSGD);

        model.swap(obj.model);
        opt.swap(obj.opt);
        mus.swap(obj.mus);
        rhos.swap(obj.rhos);
        noises.swap(obj.noises);
    }
}
