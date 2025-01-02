/*
 * Copyright 2024-2025 Weibo He.
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

#include "Vegas.h"

namespace Physica::Core {
    /**
     * \class LnVegas targets divergent integrals and computes logarithmic observations.
     */
    template<Scalar T>
    class LnVegas : public Vegas<T> {
        using This = LnVegas<T>;
        using Base = Vegas<T>;
        using typename Base::ValueType;
        using typename Base::RealValue;
        using typename Base::LossMatrix;
        using typename Base::CountArray;
    protected:
        using Base::lossMat;
        using Base::counts;

        using Base::means;
        using Base::vars;
        using Base::loss;
    public:
        using Base::Base;
        LnVegas(const This&) = default;
        LnVegas(This&&) noexcept = default;
        ~LnVegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void warmup(Functor func, int numWarm);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void integral(Functor lnFunc);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        [[nodiscard]] RealValue calcGridLoss(Functor func) const;

        [[nodiscard]] T calcLnMean() const;
        [[nodiscard]] T calcLnDevia() const;
        [[nodiscard]] T calcLnVar() const;
        [[nodiscard]] T calcSquaredChi() const;

        void swap(This& __restrict obj) noexcept { Base::swap(obj); }
        /* Getters */
        using Base::getPointGrid;
        using Base::getNumPoint;
        using Base::getDim;
        using Base::getNumRefine;
        using Base::getNumSample;
    private:
        using Base::calcMean;
        using Base::calcDevia;
        using Base::calcVar;

        template<class Functor, RandomGenerator R, class Executor>
        std::pair<T, T> trialIntegral(Functor lnFunc);
    };

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void LnVegas<T>::warmup(Functor lnFunc, int numWarm) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");
        assert(numWarm >= 0 && "[Error]: Invalid param");

        for (int _ = 0; _ < numWarm; ++_) {
            Base::pre_trial();
            trialIntegral<Functor, R, Executor>(lnFunc);
            Base::template refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void LnVegas<T>::integral(Functor lnFunc) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        for (int refine = 0; refine < getNumRefine(); ++refine) {
            Base::pre_trial();
            auto pair = trialIntegral<Functor, R, Executor>(lnFunc);
            means[refine] = std::move(pair.first);
            vars[refine] = std::move(pair.second);
            loss[refine] = Base::calcGridLossImpl();
            Base::template refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    LnVegas<T>::RealValue LnVegas<T>::calcGridLoss(Functor func) const {
        Base::pre_trial();
        trialIntegral<Functor, R, Executor>(func);
        return Base::calcGridLossImpl();
    }

    template<Scalar T>
    T LnVegas<T>::calcLnMean() const {
        return (means - vars).lnSumExp() + calcLnVar();
    }

    template<Scalar T>
    T LnVegas<T>::calcLnDevia() const {
        return RealValue(0.5) * calcLnVar();
    }

    template<Scalar T>
    T LnVegas<T>::calcLnVar() const {
        return -(-vars).lnSumExp();
    }

    template<Scalar T>
    T LnVegas<T>::calcSquaredChi() const {
        if (getNumRefine() == 1)
            return 1;
        const T lnMean = calcLnMean();
        VectorND<T> buffer = exp(means - lnMean);
        buffer = ln(abs(buffer - RealValue(1))) + lnMean;
        buffer = RealValue(2) * buffer - vars;
        return exp(buffer).sum() / RealValue(getNumRefine() - 1);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    std::pair<T, T> LnVegas<T>::trialIntegral(Functor lnFunc) {
        const int numSample = getNumSample();
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, getNumPoint() - 2);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            const auto pair = Base::template sample<R>(indexes.data_ptr(n * getDim()));
            const auto& x = pair.first;
            const auto& deltas = pair.second;
            const T lny = lnFunc(x);
            const T lnxy = lny + ln(deltas).sum() + RealValue(getDim()) * ln(RealValue(getNumPoint()));
            samples[n] = lnxy;
        }, numSample, Executor::getNumThread()).wait();

        ValueType maxSample;
        if constexpr (T::isComplex)
            maxSample = samples.reals().max().value();
        else
            maxSample = samples.max().value(); // Real LnVegas assumes f(x) > 0, so ln(f(x)) is defined
        samples = exp(samples - maxSample);

        T mean = 0, var = 0;
        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(var, mean, n, xy);

            const auto l = std::max(xy.value().squaredNorm(), RealValue(std::numeric_limits<T>::min()));
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(lossMat(index, i), counts[index][i], l);
                counts[index][i] += 1;
            }
        }
        mean = ln(mean) + maxSample;
        var = ln(var) + RealValue(2) * maxSample;
        return std::make_pair(std::move(mean), std::move(var));
    }
}
