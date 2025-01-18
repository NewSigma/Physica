/*
 * Copyright 2023-2025 Weibo He.
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
#include "Layer/LayerBase.h"

namespace Physica::Core {
    template<class Derived>
    class SeqNet : public LayerBase<Derived> {
        using This = SeqNet<Derived>;
        using Base = LayerBase<Derived>;
    public:
        using typename Base::ScalarType;
        using Base::IsTrain;
    protected:
        using Tv = ScalarType::ValueType;
    public:
        ~SeqNet() = default;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomType, class Executor>
        void train_step(const Dataset& dataset, Optimizer& opt);
        template<class Dataset, class Optimizer, class RandomType, class Executor>
        void train_step_for(int64_t numStep, const Dataset& dataset, Optimizer& opt);

        template<class Dataset>
        [[nodiscard]] auto loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;
    protected:
        SeqNet() = default;
        SeqNet(const SeqNet&) = default;
        SeqNet(SeqNet&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomType, class Executor>
    void SeqNet<Derived>::train_step(const Dataset& dataset, Optimizer& opt) {
        static_assert(IsTrain, "[Error]: train_step must be called under training mode");
        const Tv mean_grad = reciprocal(Tv(opt.getBatchSize()));
        if constexpr (std::same_as<Executor, SequentialExecutor>) {
            const auto indexes = RandomType::random_int(opt.getBatchSize(), 0, dataset.getSize() - 1);
            for (auto index : indexes)
                loss<Dataset>(dataset, index).reverse(mean_grad);
        }
        else {
            const size_t numThread = Executor::getNumThread();
            const size_t batchSizePerThread = (opt.getBatchSize() + numThread - 1) / numThread;
            std::mutex mutex{};
            Executor::parallel_for([this, mean_grad, batchSizePerThread, &dataset, &mutex](size_t) {
                auto& nn = Base::getDerived();
                Derived buffer = nn;
                const auto indexes = RandomType::random_int(batchSizePerThread, 0, dataset.getSize() - 1);
                for (auto index : indexes)
                    buffer.template loss<Dataset>(dataset, index).reverse(mean_grad);
                std::unique_lock locker(mutex);
                buffer.reverse(nn);
            }, numThread, numThread).wait();
        }
        Base::step(opt);
        Base::zero_grad();
    }

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomType, class Executor>
    void SeqNet<Derived>::train_step_for(int64_t numStep, const Dataset& dataset, Optimizer& opt) {
        for (int64_t _ = 0; _ < numStep; ++_)
            train_step<Dataset, Optimizer, RandomType, Executor>(dataset, opt);
    }

    template<class Derived>
    template<class Dataset>
    [[nodiscard]] auto SeqNet<Derived>::loss(const Dataset& dataset) const -> ScalarType {
        static_assert(!IsTrain, "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        ScalarType result = 0;
        for (size_t i = 0; i < size; ++i)
            toNextMean(result, i, loss(dataset, i));
        return result;
    }
}
