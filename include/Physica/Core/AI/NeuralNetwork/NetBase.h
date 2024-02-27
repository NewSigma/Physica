/*
 * Copyright 2023-2024 WeiBo He.
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

#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/MultiPrecision/AutoDiffGuard.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Layer/LayerBase.h"
#include "Loss.h"

namespace Physica::Core {
    template<class Derived>
    class NetBase : public LayerBase<Derived> {
        using This = NetBase<Derived>;
        using Base = LayerBase<Derived>;
        using TraitsType = typename Internal::Traits<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using typename Base::InputType;
        using typename Base::OutputType;
        using Base::IsTrainMode;
        using LossType = typename Loss<ScalarType>::LossType;
    private:
        using NetGuardType = typename std::conditional<IsTrainMode, AutoDiffGuard<ScalarType>, PlainStruct<void>>::type;

        NetGuardType net_guard;
    public:
        ~NetBase() = default;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomPoolType, class Executor>
        void train_step(const Dataset& dataset, Optimizer& opt);

        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;
        [[nodiscard]] size_t classify(const InputType& input) const;

        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        NetBase() = default;
        NetBase(const NetBase&);
        NetBase(NetBase&&) noexcept = default;
        /* Operators */
        NetBase& operator=(NetBase obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void swap(NetBase& __restrict obj) noexcept;
    };

    template<class Derived>
    NetBase<Derived>::NetBase(const NetBase&) : NetBase() {}

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomPoolType, class Executor>
    void NetBase<Derived>::train_step(const Dataset& dataset, Optimizer& opt) {
        static_assert(IsTrainMode, "[Error]: train_step must be called under training mode");
        using TracerType = typename ScalarType::TracerType;
        if constexpr (std::is_same<Executor, SequentialExecutor>::value) {
            auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
            auto& gen = RandomPoolType::getInstance().getGen();
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
                auto& gen = RandomPoolType::getInstance().getGen();
                for (size_t _ = 0; _ < batchSizePerThread; ++_) {
                    AutoDiffGuard<ScalarType> guard{};
                    const size_t index = dist(gen);
                    tempNet.template loss<Dataset>(dataset, index).reverse_to(guard.getNode());
                }
                auto& tracer = TracerType::getInstance();
                std::lock_guard locker(mutex);
                tracer.reverse(); // TODO: reverse on net_tracer should lock mutex of net_tracer
            }, numThread, numThread).wait();
        }
        opt.step();
        TracerType::getInstance().zero_grad_to(net_guard.getNode());
    }

    template<class Derived>
    template<class Dataset>
    [[nodiscard]] typename NetBase<Derived>::ScalarType NetBase<Derived>::loss(const Dataset& dataset) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        ScalarType result = 0;
        for (size_t i = 0; i < size; ++i)
            toNextMean(result, i, loss(dataset, i));
        return result;
    }

    template<class Derived>
    size_t NetBase<Derived>::classify(const InputType& input) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const OutputType output = Base::forward(input);
        PlainScalar max = output[0];
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
    inline void NetBase<Derived>::swap(NetBase& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        if constexpr (IsTrainMode)
            net_guard.swap(obj.net_guard);
    }
}
