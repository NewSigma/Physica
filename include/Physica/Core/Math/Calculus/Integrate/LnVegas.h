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
        using This = LnVegas<T>;
        using Base = Vegas<T>;
        using typename Base::ValueType;
        using typename Base::LossMatrix;

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
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        [[nodiscard]] ValueType calcGridLoss(Functor func) const;

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
        std::pair<T, T> trialIntegral(LossMatrix& losses, int refine, Functor lnFunc) const;
    };

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void LnVegas<T>::lnIntegral(Functor lnFunc) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        LossMatrix losses(getNumPoint() - 1, getDim(), 0);
        for (int refine = 0; refine < getNumRefine(); ++refine) {
            auto pair = trialIntegral<Functor, R, Executor>(losses, refine, lnFunc);
            means[refine] = std::move(pair.first);
            vars[refine] = std::move(pair.second);

            Base::template refineGrid<Executor>(losses);
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    LnVegas<T>::ValueType LnVegas<T>::calcGridLoss(Functor func) const {
        LossMatrix losses(getNumPoint() - 1, getDim(), 0);
        trialIntegral<Functor, R, Executor>(losses, 0, func);
        return Base::calcGridLossImpl(losses);
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
            return 1;
        const T lnMean = calcLnMean();
        VectorND<T> buffer = exp(means - lnMean);
        buffer = ln(abs(buffer - T(1))) + lnMean;
        buffer = T(2) * buffer - vars;
        return exp(buffer).sum() / T(getNumRefine() - 1);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    std::pair<T, T> LnVegas<T>::trialIntegral(LossMatrix& losses, int refine, Functor lnFunc) const {
        const int numSample = getNumSample();
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, losses.getRow() - 1);
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
            const T lnxy = lny + ln(deltaX).sum() + T(getDim()) * ln(T(losses.getRow()));
            samples[n] = lnxy;
        }, numSample, Executor::getNumThread()).wait();

        Array<Array<int>> counts(losses.getRow(), getDim(), 0);
        const T maxSample = samples.max();
        samples = exp(samples - maxSample);
        T mean = 0, var = 0;
        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(var, mean, n, xy);
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(losses(index, i), counts[index][i], square(xy.getValue()));
                counts[index][i] += 1;
            }
        }
        mean = ln(mean) + maxSample;
        var = ln(var) + T(2) * maxSample;
        return std::make_pair(std::move(mean), std::move(var));
    }
}
