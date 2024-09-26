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

#include "SimpleNet.h"
#include "Loss.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<SimpleNet<Derived>> : public device_obj<LayerBase<Derived>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = SimpleNet<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<Derived>>;
    public:
        using typename Base::PlainScalar;
        using typename Base::ScalarType;
        using typename Base::InputType;
        using typename Base::OutputType;
        using Base::IsTrainMode;
        using LossType = typename device_obj<Loss<ScalarType>>::LossType;
    private:
        using DiffGuardType = typename std::conditional<IsTrainMode, AutoDiffGuard<device_obj<ScalarType>>, PlainStruct<void>>::type;

        DiffGuardType net_guard;
    public:
        device_obj(const device_obj&) = delete;
        ~device_obj() = default;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomPoolType>
        void train_step(const Dataset& dataset, Optimizer& opt);

        template<class Dataset>
        [[nodiscard]] LossType loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] LossType loss(const Dataset& dataset) const;
        [[nodiscard]] size_t classify(const InputType& input) const;
        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        device_obj() = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(device_obj& __restrict obj) noexcept;
    };

    template<class Derived>
    template<class Dataset, class Optimizer, class RandomPoolType>
    void device_obj<SimpleNet<Derived>>::train_step(const Dataset& dataset, Optimizer& opt) {
        static_assert(IsTrainMode, "[Error]: train_step must be called under training mode");
        using TracerType = typename device_obj<ScalarType>::TracerType;

        auto dist = std::uniform_int_distribution<size_t>(0, dataset.getSize() - 1);
        auto& gen = RandomPoolType::getInstance().getGen();
        for (unsigned int _ = 0; _ < opt.getBatchSize(); ++_) {
            DiffGuardType guard{};
            const size_t index = dist(gen);
            loss<Dataset>(dataset, index).reverse_to(guard);
        }
        opt.step();
        TracerType::getInstance().zero_grad_to(net_guard);
    }

    template<class Derived>
    template<class Dataset>
    typename device_obj<SimpleNet<Derived>>::LossType device_obj<SimpleNet<Derived>>::loss(const Dataset& dataset) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        LossType result = 0;
        for (size_t i = 0; i < size; ++i)
            toNextMean(result, i, loss(dataset, i));
        return result;
    }

    template<class Derived>
    size_t device_obj<SimpleNet<Derived>>::classify(const InputType& input) const {
        static_assert(!IsTrainMode, "[Error]: It is suggested using eval mode to reduce memory use");
        const auto output = Base::forward(input).toHost();
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
    void device_obj<SimpleNet<Derived>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        net_guard.swap(obj.net_guard);
    }
}
