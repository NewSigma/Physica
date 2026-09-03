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

#include "Layer/LayerMixin.h"

namespace Physica {
    template<class Derived>
    class SeqNet : public LayerMixin<Derived> {
        using This = SeqNet<Derived>;
        using Base = LayerMixin<Derived>;
    protected:
        using T = Base::T;
        using Tv = Base::Tv;
    public:
        ~SeqNet() = default;
        /* Operations */
        template<RNG R, ExecutePolicy P>
        void train_step(int batchSize, const auto& dataset);
        template<RNG R, ExecutePolicy P>
        void train_step_for(int64_t numStep, int batchSize, const auto& dataset);

        [[nodiscard]] auto loss(const auto& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        [[nodiscard]] T loss(const auto& dataset) const;
    protected:
        SeqNet() = default;
        SeqNet(const SeqNet&) = default;
        SeqNet(SeqNet&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    template<RNG R, ExecutePolicy P>
    void SeqNet<Derived>::train_step(int batchSize, const auto& dataset) {
        static_assert(Base::isTraining(), "[Error]: train_step must be called under training mode");
        const Tv mean_grad = reciprocal(Tv(batchSize));
        if constexpr (P == Sequential) {
            const auto indices = R::random_int(batchSize, 0, dataset.getSize() - 1);
            for (auto index : indices)
                loss(dataset, index).reverse(mean_grad);
        }
        else {
            const int numThread = ThreadPool::getInstance().getNumThreads();
            const int batchSizePerThread = (batchSize + numThread - 1) / numThread;
            std::mutex mutex{};
            parallel_for<P>([this, mean_grad, batchSizePerThread, &dataset, &mutex](size_t) {
                auto& nn = Base::getDerived();
                Derived buffer = nn;
                const auto indices = R::random_int(batchSizePerThread, 0, dataset.getSize() - 1);
                for (auto index : indices)
                    buffer.loss(dataset, index).reverse(mean_grad);
                std::unique_lock locker(mutex);
                buffer.reverse(nn);
            }, numThread, numThread).wait();
        }
        Base::step();
        Base::zero_grad();
    }

    template<class Derived>
    template<RNG R, ExecutePolicy P>
    void SeqNet<Derived>::train_step_for(int64_t numStep, int batchSize, const auto& dataset) {
        for (int64_t _ = 0; _ < numStep; ++_)
            train_step<R, P>(batchSize, dataset);
    }

    template<class Derived>
    auto SeqNet<Derived>::loss(const auto& dataset) const -> T {
        static_assert(Base::isInfering(), "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        T result = 0;
        for (size_t i = 0; i < size; ++i)
            result.toNextMean(i, loss(dataset, i));
        return result;
    }
}
