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

#include <Physica/Core/AI/NeuralNetwork/Layer/LayerBase.h>

namespace Physica::Core {
    template<class Derived>
    class ConservedFieldNet : public LayerBase<Derived> {
        using This = ConservedFieldNet<Derived>;
        using Base = LayerBase<Derived>;
        using TraitsType = Traits<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using typename Base::InputType;
        using typename Base::OutputType;

        using device_obj_type = device_obj<LinearLayer<ScalarType, true>>;
        using DiffScalar1 = Differentiable<PlainScalar, DiffMode::Reverse, 1>;
        using LossType = typename Loss<ScalarType>::LossType;
        constexpr static bool IsTrainMode = ScalarType::Order == 2;
        static_assert(std::is_same<OutputType, ScalarType>::value, "[Error]: Output is energy, which should be a scalar");
    private:
        using DiffGuard = typename std::conditional<IsTrainMode, AutoDiffGuard<ScalarType>, PlainStruct<void>>::type;

        DiffGuard diffGuard;
    public:
        ConservedFieldNet(const ConservedFieldNet&) = delete;
        ~ConservedFieldNet() = default;
        /* Operations */
        template<class Dataset, class Optimizer, class RandomPoolType, class Executor>
        void train_step(const Dataset& dataset, Optimizer& opt);

        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        template<class Dataset>
        [[nodiscard]] ScalarType loss(const Dataset& dataset) const;

        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        ConservedFieldNet() = default;
        ConservedFieldNet(ConservedFieldNet&&) noexcept = default;
        /* Operators */
        ConservedFieldNet& operator=(ConservedFieldNet obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void swap(ConservedFieldNet& __restrict obj) noexcept;
    };

    template<class Derived>
    inline void ConservedFieldNet<Derived>::swap(ConservedFieldNet& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        if constexpr (IsTrainMode)
            diffGuard.swap(obj.diffGuard);
    }
}
