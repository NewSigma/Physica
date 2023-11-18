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

#include "Physica/Utils/Template/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    namespace Internal {
        template<class T> class Traits;
    }

    template<class Derived>
    class LayerBase : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        using PlainScalar = typename ScalarType::PlainScalar;
        using VectorType = Vector<ScalarType>;
    public:
        ~LayerBase() = default;
        /* Operations */
        [[nodiscard]] VectorType forward(const VectorType& x) const { return Base::getDerived().forward(x); }
        template<class Optimizer>
        void opt_step(const Optimizer& opt) { Base::getDerived().opt_step(opt); }
        [[nodiscard]] Derived copy() const { return Base::getDerived().copy(); }
        /* Getters */
        [[nodiscard]] constexpr static bool isTrainMode() { return Internal::Traits<ScalarType>::isDifferentiable; }
    protected:
        LayerBase() = default;
        LayerBase(const LayerBase&) = default;
        LayerBase(LayerBase&&) noexcept = default;
        /* Operators */
        LayerBase& operator=(const LayerBase&) = default;
        LayerBase& operator=(LayerBase&&) noexcept = default;
    };
}
