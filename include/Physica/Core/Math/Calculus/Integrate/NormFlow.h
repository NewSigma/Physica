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

#include "Physica/Core/ML/NeuralNetwork/Layer/LayerBase.h"
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "IntegrateImpl/AdaptiveBase.h"

namespace Physica {
    /**
     * Monte-carlo integration with neural importance sampling
     *
     * Reference:
     * [1] ACM Transactions on Graphics, 38(5) 1-19 (2019); https://doi.org/10.1145/3341156
     */
    template<Scalar T, bool TakeLn>
    class NormFlow : public AdaptiveBase<typename std::conditional<ReverseDiff<T>, typename T::ValueType, T>::type, TakeLn> {
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using FuncValue = std::conditional<ReverseDiff<T>, Tv, T>::type;
        using This = NormFlow<T, TakeLn>;
        using Base = AdaptiveBase<FuncValue, TakeLn>;
    protected:
        using Base::from;
        using Base::to;
        using Base::means;
        using Base::vars;
        using Base::loss;
    public:
        using Base::Base;
        NormFlow() = default;
        NormFlow(const This&) = default;
        NormFlow(This&&) noexcept = default;
        ~NormFlow() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net, RNG R, class Executor = SeqExecutor>
        Trv warmup(Net& nn, int numWarm);
        template<DNN Net, RNG R, class Executor = SeqExecutor>
        void integral(Net& nn);
        [[nodiscard]] T calcSquaredChi(int from = 0) const { return Base::calcSquaredChiImpl<false>(from); }

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getDim;
    private:
        VectorND<Trv> makeMask() const;
        template<DNN Net, RNG R, class Executor>
        Trv trial_normal(Net& nn, FuncValue& mean, FuncValue& var);
        template<DNN Net, RNG R, class Executor>
        Trv trial_ln(Net& nn, FuncValue& mean, FuncValue& var);

        Trv calcLoss_normal(FuncValue y, const CoDiff<T>& lnJ);
        template<class Executor>
        Trv calcLoss_ln(Trv lnVolume, VectorND<FuncValue>& samples, Array<CoDiff<T>>& lnJs, const VectorND<T>& lnJv);
    };

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    auto NormFlow<T, TakeLn>::warmup(Net& nn, int numWarm) -> Trv {
        using CallResult = std::invoke_result<Net, VectorND<Trv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<FuncValue, CallResult>, "[Error]: Inconsistent ScalarType");

        Trv result;
        for (int _ = 0; _ < numWarm; ++_) {
            FuncValue mean = 0, var = 0;
            if constexpr (TakeLn)
                result = trial_ln<Net, R, Executor>(nn, mean, var);
            else
                result = trial_normal<Net, R, Executor>(nn, mean, var);

            nn.step();
            nn.zero_grad();
        }
        return result;
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    void NormFlow<T, TakeLn>::integral(Net& nn) {
        using CallResult = std::invoke_result<Net, VectorND<Trv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<FuncValue, CallResult>, "[Error]: Inconsistent ScalarType");

        const int numRefine = Base::getNumRefine();
        FuncValue mean = 0, var = 0;
        for (int refine = 0; refine < numRefine; ++refine) {
            if constexpr (TakeLn)
                loss[refine] = trial_ln<Net, R, Executor>(nn, mean, var);
            else
                loss[refine] = trial_normal<Net, R, Executor>(nn, mean, var);
            means[refine] = mean;
            vars[refine] = var;

            nn.step();
            nn.zero_grad();
        }
    }

    template<Scalar T, bool TakeLn>
    void NormFlow<T, TakeLn>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
    }

    template<Scalar T, bool TakeLn>
    auto NormFlow<T, TakeLn>::makeMask() const -> VectorND<Trv> {
        VectorND<Trv> mask(getDim(), 0);
        mask.head(getDim() / 2) = Trv(1);
        return mask;
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    auto NormFlow<T, TakeLn>::trial_normal(Net& nn, FuncValue& mean, FuncValue& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;
        const auto volumn = (to - from).prod();
        const auto mask = makeMask();

        VectorND<Trv> x(getDim());
        Trv loss = 0;
        if constexpr (std::same_as<SeqExecutor, Executor>) {
            for (int i = 0; i < numSample; ++i) {
                x.template random_uniform<R>();
                const auto lnJ = nn.forward(x, mask);
                x = from + hadamard(coeff, x);
                FuncValue y = nn(x);
                FuncValue sample = y * exp(lnJ.value()) * volumn;
                toNextVariance(var, mean, i, sample);
                loss += calcLoss(nn(x), lnJ);
            }
        }
        else {
            const int numThread = Executor::getNumThread();
            Array<Net> nets(numThread, nn);

            VectorND<FuncValue> samples(numSample);
            VectorND<Trv> losses(numThread, 0);
            Executor::parallel_for([&, this](size_t i) {
                const int tid = Executor::getThreadID();
                auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                const auto lnJ = nets[tid].forward(x, mask);
                x = from + hadamard(coeff, x);
                FuncValue y = nn(x);
                samples[i] = y * exp(lnJ.value()) * volumn;
                losses[tid] += calcLoss(y, lnJ);
            }, numSample, numThread).wait();
            loss = losses.sum();

            for (auto& net : nets)
                nn.reverse(net);
            mean = samples.mean();
            var = samples.variance(mean);
        }
        return loss / Trv(numSample);
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    auto NormFlow<T, TakeLn>::trial_ln(Net& nn, FuncValue& mean, FuncValue& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;
        const auto lnVolume = ln(to - from).sum();
        const auto mask = makeMask();

        VectorND<FuncValue> samples(numSample);
        Array<CoDiff<T>> lnJs(numSample);
        VectorND<T> lnJv(numSample);
        Trv loss = 0;
        if constexpr (std::same_as<SeqExecutor, Executor>) {
            VectorND<Trv> x(getDim());
            for (size_t i = 0; i < numSample; ++i) {
                x.template random_uniform<R>();
                lnJs[i] = nn.forward(x, mask);
                lnJv[i] = lnJs[i].value();
                x = from + hadamard(coeff, x);
                samples[i] = nn(x);
            }
            loss = calcLoss_ln<Executor>(lnVolume, samples, lnJs, lnJv);
        }
        else {
            const int numThread = Executor::getNumThread();
            Array<Net> nets(numThread, nn);
            Executor::parallel_for([&, this](size_t i) {
                const int tid = Executor::getThreadID();
                auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                lnJs[i] = nets[tid].forward(x, mask);
                lnJv[i] = lnJs[i].value();
                x = from + hadamard(coeff, x);
                samples[i] = nn(x);
            }, numSample, numThread).wait();
            loss = calcLoss_ln<Executor>(lnVolume, samples, lnJs, lnJv);

            for (auto& net : nets)
                nn.reverse(net);
        }

        Tv maxSample;
        if constexpr (T::isComplex)
            maxSample = samples.reals().max().value();
        else
            maxSample = samples.max().value(); // Assuming f(x) > 0, so ln(f(x)) is defined
        samples = exp(samples - maxSample);

        mean = samples.mean();
        var = samples.variance(mean);
        mean = ln(mean) + maxSample;
        var = ln(var) + Trv(2) * maxSample;
        return loss / Trv(numSample);
    }

    template<Scalar T, bool TakeLn>
    auto NormFlow<T, TakeLn>::calcLoss_normal(FuncValue y, const CoDiff<T>& lnJ) -> Trv {
        CoDiff<T> l = square(y * exp(lnJ));
        if constexpr (ReverseDiff<T>)
            l.reverse(reciprocal(Trv(Base::getNumSample())));
        return l.value();
    }

    template<Scalar T, bool TakeLn>
    template<class Executor>
    auto NormFlow<T, TakeLn>::calcLoss_ln(Trv lnVolume, VectorND<FuncValue>& samples, Array<CoDiff<T>>& lnJs, const VectorND<T>& lnJv) -> Trv {
        Trv result;
        {
            auto l = (Trv(2) * (samples + lnJv)).lnSumExp();
            if constexpr (ReverseDiff<T>)
                l.reverse();
            result = l.value();
        }

        const int numThread = Executor::getNumThread();
        const size_t numSample = Base::getNumSample();
        auto futures = Executor::parallel_for([&, lnVolume](size_t i) {
            auto& lnJ = lnJs[i];
            samples[i] += lnJ.value() + lnVolume;
            if constexpr (ReverseDiff<T>)
                lnJ.reverse_final(lnJv[i].grad());
        }, numSample, numThread);
        futures.wait();
        return result;
    }
}
