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

#include <Physica/CRTPBase.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class Derived>
    class LayerBase : public CRTPBase<Derived> {
        using Base = CRTPBase<Derived>;
        using TraitsType = Traits<Derived>;
    public:
        using ScalarType = typename TraitsType::ScalarType;
        using PlainScalar = typename ScalarType::PlainScalar;
        using InputType = typename TraitsType::InputType;
        using OutputType = typename TraitsType::OutputType;
        constexpr static bool IsTrainMode = TraitsType::IsTrainMode;
    public:
        ~LayerBase() = default;
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const { return Base::getDerived().forward(x); }
        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
    protected:
        LayerBase() = default;
        LayerBase(const LayerBase&) = default;
        LayerBase(LayerBase&&) noexcept = default;
        /* Operators */
        LayerBase& operator=(const LayerBase&) = default;
        LayerBase& operator=(LayerBase&&) noexcept = default;
    };
}
