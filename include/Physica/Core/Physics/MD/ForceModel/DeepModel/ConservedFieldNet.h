/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/ML/NeuralNetwork/Layer/LayerBase.h"

namespace Physica {
    template<class Derived>
    class ConservedFieldNet : public LayerBase<Derived> {
        using This = ConservedFieldNet<Derived>;
        using Base = LayerBase<Derived>;
        using TraitsType = Traits<Derived>;
    public:
        using typename Base::ScalarType;
        using Tv = ScalarType::ValueType;

        using device_obj_type = device_obj<LinearLayer<ScalarType, true>>;
        using DiffScalar1 = Diff<Tv, DiffMode::Reverse, 1>;
        using LossType = Loss<ScalarType>::LossType;
    public:
        ConservedFieldNet(const This&) = delete;
        ~ConservedFieldNet() = default;
        /* Operations */
        template<class RandomSource, ExecutePolicy P>
        void train_step(const auto& dataset, auto& optimizer) { noImpl(); }

        [[nodiscard]] ScalarType loss(const auto& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        [[nodiscard]] ScalarType loss(const auto& dataset) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ consteval static bool isTraining() noexcept { return ScalarType::Order == 2; }
    protected:
        ConservedFieldNet() = default;
        ConservedFieldNet(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
    };

    template<class Derived>
    void ConservedFieldNet<Derived>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        if constexpr (isTraining())
            diffGuard.swap(obj.diffGuard);
    }
}
