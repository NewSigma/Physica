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

#include "LayerBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<LayerBase<Derived>> : public CRTPBase<device_obj<LayerBase<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = LayerBase<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
        using device_obj_type = device_obj<Derived>;
        using TraitsType = Traits<device_obj_type>;
    public:
        using ScalarType = TraitsType::ScalarType;
        using ValueType = ScalarType::ValueType;
        using OutputType = TraitsType::OutputType;
        constexpr static bool IsTrain = ScalarType::isDiffable;
        static_assert(!is_device_obj<ScalarType>::value, "[Error]: Nested device_obj<> is not allowed");
    public:
        ~device_obj() = default;
        /* Operations */
        template<class T>
        [[nodiscard]] OutputType forward(const T& x) const { return Base::getDerived().template forward<T>(x); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<LayerBase<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}
