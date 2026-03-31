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

#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "IntegrateImpl/AdaptiveBase.h"

namespace Physica {
    class CUDAExecutor;
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
        using IntegT = std::conditional<ReverseDiff<T>, Tv, T>::type;
        using TrainT = std::conditional<ReverseDiff<T>, T, Tv>::type;
        using This = NormFlow<T, TakeLn>;
        using Base = AdaptiveBase<IntegT, TakeLn>;

        template<Scalar U>
        using MatrixND = LinearLayer<T>::template MatrixND<U>;
    protected:
        using Base::from;
        using Base::to;
        using Base::means;
        using Base::vars;
        using Base::loss;
    private:
        int batchsize;
        Trv decay;
    public:
        NormFlow() = default;
        NormFlow(VectorND<Trv> from, VectorND<Trv> to, int numRefine, int numSample, int batchsize_, Trv decay_ = 0);
        NormFlow(const This&) = default;
        NormFlow(This&&) noexcept = default;
        ~NormFlow() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net, RNG R, ExecutePolicy P = Sequential>
        Trv warmup(Net& nn, int numWarm);
        template<DNN Net, RNG R, ExecutePolicy P = Sequential>
        void integral(Net& nn);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getDim;
    private:
        template<DNN Net, RNG R, ExecutePolicy P>
        Trv trial_normal(Net& nn, IntegT& mean, IntegT& var);
        template<DNN Net, RNG R, ExecutePolicy P>
        Trv trial_ln(Net& nn, IntegT& mean, IntegT& var);

        Trv calcLoss_normal(IntegT y, const CoDiff<TrainT>& lnJ);
        template<ExecutePolicy P>
        Trv calcLoss_ln(const Vector auto& samples, const VectorND<TrainT>& lnJv);
    };

    template<Scalar T, bool TakeLn>
    NormFlow<T, TakeLn>::NormFlow(VectorND<Trv> from, VectorND<Trv> to, int numRefine, int numSample, int batchsize_, Trv decay_)
            : Base(std::move(from), std::move(to), numRefine, (numSample + batchsize_ - 1) / batchsize_ * batchsize_)
            , batchsize(batchsize_)
            , decay(decay_) {
        assert(0 < batchsize && batchsize <= numSample && "[Error]: Invalid batchsize and cannot auto fix");
        assert(!decay.isNegative());
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, ExecutePolicy P>
    auto NormFlow<T, TakeLn>::warmup(Net& nn, int numWarm) -> Trv {
        using CallResult = std::invoke_result<Net, VectorND<Trv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<IntegT, CallResult>, "[Error]: Inconsistent ScalarType");

        Trv result;
        for (int _ = 0; _ < numWarm; ++_) {
            IntegT mean = 0, var = 0;
            if constexpr (TakeLn)
                result = trial_ln<Net, R, P>(nn, mean, var);
            else
                result = trial_normal<Net, R, P>(nn, mean, var);

            nn.step();
            nn.zero_grad();
        }
        return result;
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, ExecutePolicy P>
    void NormFlow<T, TakeLn>::integral(Net& nn) {
        using CallResult = std::invoke_result<Net, VectorND<Trv>>::type;
        static_assert(Scalar<CallResult>, "[Error]: Integrand should embed into network");
        static_assert(std::same_as<IntegT, CallResult>, "[Error]: Inconsistent ScalarType");

        const int numRefine = Base::getNumRefine();
        IntegT mean = 0, var = 0;
        for (int refine = 0; refine < numRefine; ++refine) {
            if constexpr (TakeLn)
                loss[refine] = trial_ln<Net, R, P>(nn, mean, var);
            else
                loss[refine] = trial_normal<Net, R, P>(nn, mean, var);
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
        std::swap(batchsize, obj.batchsize);
        decay.swap(obj.decay);
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, ExecutePolicy P>
    auto NormFlow<T, TakeLn>::trial_normal(Net& nn, IntegT& mean, IntegT& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;

        VectorND<Trv> x(getDim());
        Trv loss = 0;
        if constexpr (P == Sequential) {
            for (int i = 0; i < numSample; ++i) {
                x.template random_uniform<R>();
                const auto lnJ = nn.forward(x);
                x = from + hadamard(coeff, x);
                IntegT y = nn(x);
                IntegT sample = y * exp(lnJ.value());
                var.toNextVariance(mean, i, sample);
                loss += calcLoss_normal(nn(x), lnJ);
            }
            loss /= Trv(numSample);
        }
        else {
            static_assert(P == Thread, "[Error]: Unsupported executor");
            const int numThread = ThreadPool::getInstance().getNumThreads();
            Array<Net> nets(numThread, nn);

            VectorND<IntegT> samples(numSample);
            VectorND<Trv> losses(numThread, 0);
            parallel_for<P>([&, this](size_t i) {
                // TODO: Use a chained wait to avoid tid use
                const int tid = ThreadPool::getInstance().getThreadID();
                auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                const auto lnJ = nets[tid].forward(x);
                x = from + hadamard(coeff, x);
                IntegT y = nn(x);
                samples[i] = y * exp(lnJ.value());
                losses[tid] += calcLoss_normal(y, lnJ);
            }, numSample, numThread).wait();

            for (auto& net : nets)
                nn.reverse(net);
            loss = losses.mean();
            mean = samples.mean();
            var = samples.variance();
        }
        const auto volume = (to - from).prod();
        mean *= volume;
        var *= square(volume);
        return loss - Trv(1);
    }

    template<Scalar T, bool TakeLn>
    template<DNN Net, RNG R, ExecutePolicy P>
    auto NormFlow<T, TakeLn>::trial_ln(Net& nn, IntegT& mean, IntegT& var) -> Trv {
        const int numSample = Base::getNumSample();
        const VectorND<Trv> coeff = to - from;
        const auto lnVolume = ln(coeff).sum();

        VectorND<IntegT> samples(numSample);
        Trv loss = 0;
        if constexpr (P == GPU) {
            const auto from_d = from.toDeviceAsync();
            const auto coeff_d = coeff.toDeviceAsync();
            const int numBatch = numSample / batchsize;

            VectorND<TrainT> lnJv(batchsize);
            device_obj<MatrixND<Tv>> x(getDim(), batchsize);
            for (int i = 0; i < numBatch; ++i) {
                x.template random_uniform<R>();
                const auto lnJs = nn.forward(x);
                auto seg = samples.segment(i * batchsize, (i + 1) * batchsize);

                lnJs.toHostAsync(lnJv);
                x = from_d + hadamard(coeff_d, x);
                constexpr bool DeviceEval = DeviceObj<decltype(nn(x))>;
                if constexpr (DeviceEval)
                    nn(x).toHost(seg);
                else
                    seg = nn(x);

                loss = calcLoss_ln<P>(seg, lnJv);
                if constexpr (ReverseDiff<T>)
                    lnJs.reverse(lnJv.grads().toDeviceAsync());
                seg += lnJv.values() + lnVolume;
            }
        }
        else {
            Array<CoDiff<T>> lnJs(numSample);
            VectorND<T> lnJv(batchsize);
            Array<Net> nets;
            if constexpr (P == Sequential) {
                VectorND<Trv> x(getDim());
                for (int i = 0; i < numSample; ++i) {
                    x.template random_uniform<R>();
                    lnJs[i] = nn.forward(x);
                    lnJv[i] = lnJs[i].value();
                    x = from + hadamard(coeff, x);
                    samples[i] = nn(x);
                }
            }
            else {
                static_assert(P == Thread, "[Error]: Unsupported executor");
                const int numThread = ThreadPool::getInstance().getNumThreads();
                nets.resize(numThread, nn);
                parallel_for<P>([&, this](size_t i) {
                    // TODO: Use a chained wait to avoid tid use
                    const int tid = ThreadPool::getInstance().getThreadID();
                    auto x = VectorND<Trv>::template random_uniform<R>(getDim());
                    lnJs[i] = nets[tid].forward(x);
                    lnJv[i] = lnJs[i].value();
                    x = from + hadamard(coeff, x);
                    samples[i] = nn(x);
                }, numSample, numThread).wait();

                for (auto& net : nets)
                    nn.reverse(net);
            }

            loss = calcLoss_ln<P>(samples, lnJv);
            parallel_for<P>([&, lnVolume](size_t i) {
                auto& lnJ = lnJs[i];
                samples[i] += lnJ.value() + lnVolume;
                if constexpr (ReverseDiff<T>)
                    lnJ.reverse_final(lnJv[i].grad());
            }, numSample, 0).wait();
        }

        Tv maxSample;
        if constexpr (T::isComplex())
            maxSample = samples.reals().max().value();
        else
            maxSample = samples.max().value(); // Assuming f(x) > 0, so ln(f(x)) is defined
        samples = exp(samples - maxSample);

        mean = ln(samples.mean()) + maxSample;
        var = ln(samples.variance() + Trv(std::numeric_limits<T>::min())) + Trv(2) * maxSample;
        return loss - ln(Trv(numSample));
    }

    template<Scalar T, bool TakeLn>
    auto NormFlow<T, TakeLn>::calcLoss_normal(IntegT y, const CoDiff<TrainT>& lnJ) -> Trv {
        CoDiff<T> l = square(y * exp(lnJ)) + exp(lnJ * Trv(2)) * decay;
        if constexpr (ReverseDiff<T>)
            l.reverse(reciprocal(Trv(Base::getNumSample())));
        return l.value();
    }

    template<Scalar T, bool TakeLn>
    template<ExecutePolicy P>
    auto NormFlow<T, TakeLn>::calcLoss_ln(const Vector auto& samples, const VectorND<TrainT>& lnJv) -> Trv {
        const size_t size = samples.getLength();
        const auto mean = (samples + lnJv.values()).lnSumExp() - ln(Trv(size));
        auto l = (Trv(2) * (samples + lnJv - mean)).lnSumExp() + lnJv.squaredNorm() * (decay / Trv(size));
        if constexpr (ReverseDiff<T>)
            l.reverse();
        return l.value();
    }
}
