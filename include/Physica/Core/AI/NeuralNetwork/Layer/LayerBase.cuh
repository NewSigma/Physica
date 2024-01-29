/*
 * Copyright 2024 WeiBo He.
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

#include "LayerBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LayerBase<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using host_obj = LayerBase<Derived>;
        using This = device_obj<host_obj>;
        using Base = Utils::CRTPBase<device_obj<Derived>>;
    public:
        constexpr static bool IsTrainMode = host_obj::IsTrainMode;
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        using PlainScalar = typename ScalarType::PlainScalar;
        using VectorType = Vector<ScalarType>;
        using DiffVector = Differentiable<Vector<PlainScalar>, DiffMode::Reverse>;
        using InputType = device_obj<typename std::conditional<IsTrainMode, DiffVector, VectorType>::type>;
        using OutputType = InputType;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const { return Base::getDerived().forward(x); }
        [[nodiscard]] device_obj<Derived> copy() const { return Base::getDerived().copy(); }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
    };
}
