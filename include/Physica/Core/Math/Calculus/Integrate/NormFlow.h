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

#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "Physica/Core/ML/NeuralNetwork/Layer/LayerBase.h"
#include "IntegrateImpl/AdaptiveBase.h"

namespace Physica {
    /**
     * Monte-carlo integration with neural importance sampling
     *
     * Reference:
     * [1] ACM Transactions on Graphics, 38(5) 1-19 (2019); https://doi.org/10.1145/3341156
     */
    template<Scalar T, bool TakeLn>
    class NormFlow : public AdaptiveBase<typename T::ValueType, TakeLn> {
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using This = NormFlow<T, TakeLn>;
        using Base = AdaptiveBase<Tv, TakeLn>;
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
        Trv trial_normal(Net& nn, Tv& mean, Tv& var);
        template<DNN Net, RNG R, class Executor>
        Trv trial_ln(Net& nn, Tv& mean, Tv& var);

        Trv calcLoss(Trv y, const CoDiff<T>& lnJ);
    };

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    auto NormFlow<T, TakeLn>::warmup(Net& nn, int numWarm) -> Trv {
        using CallResult = std::invoke_result<Net, VectorND<Tv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<T, typename Traits<Net>::ScalarType>, "[Error]: Inconsistent ScalarType");

        T mean = 0, var = 0;
        Trv result;
        for (int _ = 0; _ < numWarm; ++_) {
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
        using CallResult = std::invoke_result<Net, VectorND<Tv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<T, typename Traits<Net>::ScalarType>, "[Error]: Inconsistent ScalarType");

        const int numRefine = Base::getNumRefine();
        for (int refine = 0; refine < numRefine; ++refine) {
            if constexpr (TakeLn)
                loss[refine] = trial_ln<Net, R, Executor>(nn, means[refine], vars[refine]);
            else
                loss[refine] = trial_normal<Net, R, Executor>(nn, means[refine], vars[refine]);

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
    auto NormFlow<T, TakeLn>::trial_normal(Net& nn, Tv& mean, Tv& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;
        const Trv volumn = (to - from).prod();
        const auto mask = makeMask();

        VectorND<Trv> x(getDim());
        Trv loss = 0;
        if constexpr (std::same_as<SeqExecutor, Executor>) {
            for (int i = 0; i < numSample; ++i) {
                x.template random_uniform<R>();
                const auto lnJ = nn.forward(x, mask);
                x = from + hadamard(coeff, x);
                Trv y = nn(x);
                Trv sample = y * exp(lnJ.value()) * volumn;
                toNextVariance(var, mean, i, sample);
                loss += calcLoss(nn(x), lnJ);
            }
        }
        else {
            const int numThread = Executor::getNumThread();
            Array<Net> nets(numThread, nn);

            VectorND<Trv> samples(numSample);
            Executor::parallel_for([&, this](size_t i) {
                const int tid = Executor::getThreadID();
                auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                const auto lnJ = nets[tid].forward(x, mask);
                x = from + hadamard(coeff, x);
                Trv y = nn(x);
                samples[i] = y * exp(lnJ.value()) * volumn;
                loss += calcLoss(y, lnJ);
            }, numSample, numThread).wait();

            for (auto& net : nets)
                nn.reverse(net);
            mean = Physica::mean(samples);
            var = Physica::variance(samples, mean);
        }
        return loss / Trv(numSample);
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, class Executor>
    auto NormFlow<T, TakeLn>::trial_ln(Net& nn, Tv& mean, Tv& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;
        const Trv lnVolumn = ln(to - from).sum();
        const auto mask = makeMask();

        VectorND<Trv> samples(numSample);
        Trv loss = 0;
        if constexpr (std::same_as<SeqExecutor, Executor>) {
            VectorND<Trv> x(getDim());
            for (size_t i = 0; i < numSample; ++i) {
                x.template random_uniform<R>();
                const auto lnJ = nn.forward(x, mask);
                x = from + hadamard(coeff, x);
                Trv lny = nn(x);
                samples[i] = lny + lnJ.value() + lnVolumn;
                loss += calcLoss(lny, lnJ);
            }
        }
        else {
            const int numThread = Executor::getNumThread();
            Array<Net> nets(numThread, nn);
            Executor::parallel_for([&, this](size_t i) {
                const int tid = Executor::getThreadID();
                auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                const auto lnJ = nets[tid].forward(x, mask);
                x = from + hadamard(coeff, x);
                Trv lny = nn(x);
                samples[i] = lny + lnJ.value() + lnVolumn;
                loss += calcLoss(lny, lnJ);
            }, numSample, numThread).wait();

            for (auto& net : nets)
                nn.reverse(net);
        }

        Tv maxSample;
        if constexpr (T::isComplex)
            maxSample = samples.reals().max().value();
        else
            maxSample = samples.max().value(); // Real LnVegas assumes f(x) > 0, so ln(f(x)) is defined
        samples = exp(samples - maxSample);

        mean = Physica::mean(samples);
        var = Physica::variance(samples, mean);
        mean = ln(mean) + maxSample;
        var = ln(var) + Trv(2) * maxSample;
        return loss / Trv(numSample);
    }

    template<Scalar T, bool TakeLn>
    auto NormFlow<T, TakeLn>::calcLoss(Trv y, const CoDiff<T>& lnJ) -> Trv {
        CoDiff<T> l;
        if constexpr (TakeLn) {
            Trv lny = y;
            l = exp(Tv(2) * (lny + lnJ));
        }
        else
            l = square(y * exp(lnJ));

        if constexpr (ReverseDiff<T>)
            l.reverse(reciprocal(Trv(Base::getNumSample())));
        return l.value();
    }
}
