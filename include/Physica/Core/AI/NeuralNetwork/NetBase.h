/*
 * Copyright 2023 WeiBo He.
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

#include "Layer/LayerBase.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/MultiPrecision/AutoDiffGuard.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

namespace Physica::Core {
    template<class Derived>
    class NetBase : public LayerBase<Derived> {
        using This = NetBase<Derived>;
        using Base = LayerBase<Derived>;
        using DerivedTraits = typename Internal::Traits<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using typename Base::VectorType;
        using SampleType = typename DerivedTraits::SampleType;
        using LabelType = typename DerivedTraits::LabelType;
    private:
        std::unique_ptr<AutoDiffGuard<PlainScalar>> net_guard;
    public:
        ~NetBase() = default;
        /* Operators */
        NetBase& operator=(const NetBase&) = delete;
        NetBase& operator=(NetBase&&) noexcept = delete;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomPoolType, class Executor>
        void train_step(const Dataset& dataset, const Optimizer& opt);
        template<class Optimizer>
        void opt_step(const Optimizer& opt) { Base::getDerived().opt_step(opt); }

        [[nodiscard]] ScalarType loss(const SampleType& sample, const LabelType& label) const { return Base::getDerived().loss(sample, label); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;
        [[nodiscard]] size_t classify(const VectorType& input) const;

        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        NetBase();
        NetBase(const NetBase&);
        NetBase(NetBase&&) noexcept = default;
        /* Operations */
        inline void swap(NetBase& obj) noexcept;
    };

    template<class Derived>
    NetBase<Derived>::NetBase() {
        if constexpr (Base::isTrainMode())
            net_guard = std::make_unique<AutoDiffGuard<PlainScalar>>();
    }

    template<class Derived>
    NetBase<Derived>::NetBase(const NetBase&) : NetBase() {}

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomPoolType, class Executor>
    void NetBase<Derived>::train_step(const Dataset& dataset, const Optimizer& opt) {
        constexpr bool isTrainMode = Base::isTrainMode();
        static_assert(isTrainMode, "[Error]: train_step must be called under training mode");
        if constexpr (std::is_same<Executor, SequentialExecutor>::value) {
            auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
            auto& gen = RandomPoolType::getGen();
            const auto& samples = dataset.getSamples();
            const auto& labels = dataset.getLabels();
            for (unsigned int _ = 0; _ < opt.getBatchSize(); ++_) {
                AutoDiffGuard<PlainScalar> guard{};
                const size_t index = dist(gen);
                loss(samples[index], labels[index]).reverse();
            }
        }
        else {
            const size_t batchSizePerThread = opt.getBatchSize() / Executor::getNumThread();
            auto& net_tracer = DiffTracer<PlainScalar>::getInstance();
            std::mutex mutex{};
            Executor::parallel_for([this, batchSizePerThread, &mutex, &net_tracer, &dataset](unsigned int) {
                Derived tempNet = copy();
                auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
                auto& gen = RandomPoolType::getGen();
                const auto& samples = dataset.getSamples();
                const auto& labels = dataset.getLabels();
                for (size_t _ = 0; _ < batchSizePerThread; ++_) {
                    AutoDiffGuard<PlainScalar> guard{};
                    const size_t index = dist(gen);
                    tempNet.loss(samples[index], labels[index]).reverse(guard.getTraceIndex());
                }
                auto& tracer = DiffTracer<PlainScalar>::getInstance();
                std::lock_guard locker(mutex);
                tracer.reverse(net_tracer, tracer.getNumRecord() - 1); // TODO: reverse on net_tracer should lock mutex of net_tracer
                //tracer.release(); // TODO: When should we release DiffTracer
            }, Executor::getNumThread(), Executor::getNumThread()).wait();
        }
        opt_step(opt);
        DiffTracer<PlainScalar>::getInstance().zero_grad(net_guard->getTraceIndex());
    }

    template<class Derived>
    template<class Dataset>
    [[nodiscard]] typename NetBase<Derived>::ScalarType NetBase<Derived>::loss(const Dataset& dataset) const {
        const size_t numSample = dataset.getSize();
        const auto& samples = dataset.getSamples();
        const auto& labels = dataset.getLabels();
        ScalarType result = 0;
        for (size_t i = 0; i < numSample; ++i)
            toNextMean(result, i, loss(samples[i], labels[i]));
        return result;
    }

    template<class Derived>
    size_t NetBase<Derived>::classify(const VectorType& input) const {
        const VectorType prob = softmax(forward(input));
        PlainScalar max = 0;
        size_t index = 0;
        for (size_t i = 0; i < prob.getLength(); ++i) {
            if (prob[i] > max) {
                index = i;
                max = prob[i];
            }
        }
        return index;
    }

    template<class Derived>
    inline void NetBase<Derived>::swap(NetBase& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        net_guard.swap(obj.net_guard);
    }
}
