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

#include "SeqNet.h"

namespace Physica {
    template<class Derived>
    class device_obj<SeqNet<Derived>> : public device_obj<LayerBase<Derived>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = SeqNet<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::IsTrain;
    protected:
        using Tv = ScalarType::ValueType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<class Dataset, RNG R, class Executor>
        void train_step(int batchSize, const Dataset& dataset);
        template<class Dataset, RNG R, class Executor>
        void train_step_for(int64_t numStep, int batchSize, const Dataset& dataset);

        template<class Dataset>
        [[nodiscard]] auto loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    template<class Dataset, RNG R, class Executor>
    void device_obj<SeqNet<Derived>>::train_step(int batchSize, const Dataset& dataset) {
        static_assert(IsTrain, "[Error]: train_step must be called under training mode");
        const auto indices = R::random_int(batchSize, 0, dataset.getSize() - 1);
        const Tv mean_grad = reciprocal(Tv(batchSize));
        for (auto index : indices)
            loss<Dataset>(dataset, index).reverse(mean_grad);
        Base::step();
        Base::zero_grad();
    }

    template<class Derived>
    template<class Dataset, RNG R, class Executor>
    void device_obj<SeqNet<Derived>>::train_step_for(int64_t numStep, int batchSize, const Dataset& dataset) {
        for (int64_t _ = 0; _ < numStep; ++_)
            train_step<Dataset, R, Executor>(batchSize, dataset);
    }

    template<class Derived>
    template<class Dataset>
    auto device_obj<SeqNet<Derived>>::loss(const Dataset& dataset) const -> ScalarType {
        static_assert(!IsTrain, "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        ScalarType result = 0;
        for (size_t i = 0; i < size; ++i)
            toNextMean(result, i, loss(dataset, i));
        return result;
    }
}
