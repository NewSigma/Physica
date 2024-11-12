/*
 * Copyright 2023-2024 Weibo He.
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

#include <Physica/Core/Math/Statistics/NumCharacter.h>
#include <Physica/Core/MultiPrecision/AutoDiffGuard.h>
#include <Physica/Core/Parallel/Executor/ThreadExecutor.h>
#include "Layer/LayerBase.h"
#include "Loss.h"

namespace Physica::Core {
    template<class Derived>
    class SimpleNet : public LayerBase<Derived> {
        using This = SimpleNet<Derived>;
        using Base = LayerBase<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using typename Base::InputType;
        using typename Base::OutputType;
        using Base::IsTrainMode;
        using LossType = typename Loss<ScalarType>::LossType;
    private:
        using DiffGuard = typename std::conditional<IsTrainMode, AutoDiffGuard<ScalarType>, PlainStruct<void>>::type;

        DiffGuard diffGuard;
    public:
        SimpleNet(const SimpleNet&) = delete;
        ~SimpleNet() = default;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomType, class Executor>
        void train_step(const Dataset& dataset, Optimizer& opt);

        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;
        [[nodiscard]] size_t classify(const InputType& input) const;

        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        SimpleNet() = default;
        SimpleNet(SimpleNet&&) noexcept = default;
        /* Operators */
        SimpleNet& operator=(SimpleNet obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void swap(SimpleNet& __restrict obj) noexcept;
    };

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomType, class Executor>
    void SimpleNet<Derived>::train_step(const Dataset& dataset, Optimizer& opt) {
        static_assert(IsTrainMode, "[Error]: train_step must be called under training mode");
        using TracerType = typename ScalarType::TracerType;
        if constexpr (std::is_same<Executor, SequentialExecutor>::value) {
            auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
            auto& gen = RandomType::getInstance().getGen();
            for (unsigned int _ = 0; _ < opt.getBatchSize(); ++_) {
                AutoDiffGuard<ScalarType> guard{};
                const size_t index = dist(gen);
                loss<Dataset>(dataset, index).reverse();
            }
        }
        else {
            const size_t numThread = Executor::getNumThread();
            const size_t batchSizePerThread = (opt.getBatchSize() + numThread - 1) / numThread;
            std::mutex mutex{};
            Executor::parallel_for([this, batchSizePerThread, &mutex, &dataset](unsigned int) {
                Derived tempNet = copy();
                auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
                auto& gen = RandomType::getInstance().getGen();
                for (size_t _ = 0; _ < batchSizePerThread; ++_) {
                    AutoDiffGuard<ScalarType> guard{};
                    const size_t index = dist(gen);
                    tempNet.template loss<Dataset>(dataset, index).reverse_to(guard);
                }
                auto& tracer = TracerType::getInstance();
                std::lock_guard locker(mutex);
                tracer.reverse();
            }, numThread, numThread).wait();
        }
        opt.step();
        TracerType::getInstance().zero_grad_to(diffGuard);
    }

    template<class Derived>
    template<class Dataset>
    [[nodiscard]] typename SimpleNet<Derived>::ScalarType SimpleNet<Derived>::loss(const Dataset& dataset) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        ScalarType result = 0;
        for (size_t i = 0; i < size; ++i)
            toNextMean(result, i, loss(dataset, i));
        return result;
    }

    template<class Derived>
    size_t SimpleNet<Derived>::classify(const InputType& input) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const OutputType output = Base::forward(input);
        ValueType max = output[0];
        size_t index = 0;
        for (size_t i = 1; i < output.getLength(); ++i) {
            if (output[i] > max) {
                index = i;
                max = output[i].getValue();
            }
        }
        return index;
    }

    template<class Derived>
    inline void SimpleNet<Derived>::swap(SimpleNet& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        diffGuard.swap(obj.diffGuard);
    }
}
