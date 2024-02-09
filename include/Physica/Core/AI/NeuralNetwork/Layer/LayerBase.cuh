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
    class device_obj<LayerBase<Derived>> : public Utils::CRTPBase<typename Internal::Traits<Derived>::device_obj_type> {
        using host_obj = LayerBase<Derived>;
        using This = device_obj<host_obj>;
        using device_obj_type = typename Internal::Traits<Derived>::device_obj_type;
        using Base = Utils::CRTPBase<device_obj_type>;
        using TraitsType = Internal::Traits<device_obj_type>;
    public:
        using PlainScalar = typename TraitsType::PlainScalar;
        using ScalarType = typename TraitsType::ScalarType;
        using InputType = typename TraitsType::InputType;
        using OutputType = typename TraitsType::OutputType;
        constexpr static bool IsTrainMode = TraitsType::IsTrainMode;
        static_assert(!Utils::is_device_obj<ScalarType>::value, "[Error]: Nested device_obj<> is not allowed");
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const { return Base::getDerived().forward(x); }
        [[nodiscard]] device_obj_type copy() const { return Base::getDerived().copy(); }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
    };
}
