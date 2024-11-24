/*
 * Copyright 2024 Weibo He.
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
    class LnVegas : private Vegas<T> {
        using Base = Vegas<T>;
        using This = LnVegas<T>;

        using Base::means;
        using Base::vars;
    public:
        using Base::Base;
        LnVegas(const This&) = default;
        LnVegas(This&&) noexcept = default;
        ~LnVegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void lnIntegral(Functor lnFunc);
        template<class Functor, RandomGenerator R>
        [[nodiscard]] T accessMerit(Functor func);

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
        /* Setters */
        using Base::setNumRefine;
    private:
        template<class Functor, RandomGenerator R, class Executor>
        void trialIntegral(DenseMatrix<T>& varsDevia, int refine, Functor lnFunc);
    };

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void LnVegas<T>::lnIntegral(Functor lnFunc) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        DenseMatrix<T> varsDevia(getNumPoint() - 1, getDim(), 0);
        for (int refine = 0; refine < getNumRefine(); ++refine) {
            trialIntegral<Functor, R, Executor>(varsDevia, refine, lnFunc);
            Base::template refineGrid<Executor>(varsDevia);
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R>
    T LnVegas<T>::accessMerit(Functor func) {
        DenseMatrix<T> varsDevia(getNumPoint() - 1, getDim(), 0);
        trialIntegral<Functor, R, SequentialExecutor>(varsDevia, 0, func);
        return Base::accessMeritImpl(varsDevia);
    }

    template<Scalar T>
    T LnVegas<T>::calcLnMean() const {
        return (means - vars).lnSumExp() + calcLnVar();
    }

    template<Scalar T>
    T LnVegas<T>::calcLnDevia() const {
        return T(0.5) * calcLnVar();
    }

    template<Scalar T>
    T LnVegas<T>::calcLnVar() const {
        return -(-vars).lnSumExp();
    }

    template<Scalar T>
    T LnVegas<T>::calcSquaredChi() const {
        if (getNumRefine() == 1)
            return 0;
        const T lnMean = calcLnMean();
        VectorND<T> buffer = exp(means - lnMean);
        buffer = ln(abs(buffer - T(1))) + lnMean;
        buffer = T(2) * buffer - vars;
        return exp(buffer).sum() / T(getNumRefine() - 1);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void LnVegas<T>::trialIntegral(DenseMatrix<T>& varsDevia, int refine, Functor lnFunc) {
        const int numSample = getNumSample();
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, varsDevia.getRow() - 1);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            VectorND<T> fromX(getDim());
            VectorND<T> deltaX(getDim());
            for (size_t i = 0; i < getDim(); ++i) {
                const auto& pointGrid = getPointGrid();
                const auto index = indexes[n * getDim() + i];
                fromX[i] = pointGrid(index, i);
                deltaX[i] = pointGrid(index + 1, i);
            }
            deltaX -= fromX;

            const VectorND<T> x = fromX + hadamard(deltaX, VectorND<T>::template random_uniform<R>(getDim()));
            const T lny = lnFunc(x);
            const T lnxy = lny + ln(deltaX).sum() + T(getDim()) * ln(T(varsDevia.getRow()));
            samples[n] = lnxy;
        }, numSample, Executor::getNumThread()).wait();

        Array<Array<int>> counts(varsDevia.getRow(), getDim(), 0);
        const T maxSample = samples.max();
        samples = exp(samples - maxSample);
        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(vars[refine], means[refine], n, xy);
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(varsDevia(index, i), counts[index][i], square(xy));
                counts[index][i] += 1;
            }
        }
        means[refine] = ln(means[refine]) + maxSample;
        vars[refine] = ln(vars[refine]) + T(2) * maxSample;
    }
}
