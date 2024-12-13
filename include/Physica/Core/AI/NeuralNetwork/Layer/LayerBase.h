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

#include "Physica/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica::Core {
    template<class Derived>
    class LayerBase : public CRTPBase<LayerBase<Derived>> {
        using This = LayerBase<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;
    public:
        using ScalarType = TraitsType::ScalarType;
        using ValueType = ScalarType::ValueType;
        using InputType = TraitsType::InputType;
        using OutputType = TraitsType::OutputType;
        constexpr static bool IsTrainMode = ScalarType::isDiffable;
    public:
        ~LayerBase() = default;
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const { return Base::getDerived().forward(x); }
        /* Getters */
        [[nodiscard]] size_t getInputDim() const noexcept { return Base::getDerived().getInputDim(); }
        [[nodiscard]] size_t getOutputDim() const noexcept { return Base::getDerived().getOutputDim(); }
    protected:
        LayerBase() = default;
        LayerBase(const This&) = default;
        LayerBase(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T>
    class Traits<LayerBase<T>> {
    public:
        using Derived = T;
    };
}
