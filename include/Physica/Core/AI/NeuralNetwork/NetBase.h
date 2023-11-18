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

namespace Physica::Core {
    template<class Derived>
    class NetBase : public LayerBase<Derived> {
        using Base = LayerBase<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using typename Base::VectorType;
    public:
        ~NetBase() = default;
        /* Operations */
        void init() { Base::getDerived().init(); }
        template<class Dataset, class Optimizer, class Executor>
        void train_step(const Dataset& dataset, const Optimizer& opt) { Base::getDerived().template train_step<Dataset, Optimizer, Executor>(dataset, opt); }
        [[nodiscard]] ScalarType loss(const VectorType& input, const VectorType& answer) { return Base::getDerived().loss(input, answer); }

        [[nodiscard]] size_t classify(const VectorType& input) const;
    protected:
        NetBase() = default;
        NetBase(const NetBase&) = default;
        NetBase(NetBase&&) noexcept = default;
        /* Operators */
        NetBase& operator=(const NetBase&) = default;
        NetBase& operator=(NetBase&&) noexcept = default;
    };

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
}
