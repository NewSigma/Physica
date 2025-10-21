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
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<LayerBase<Derived>> : public CRTPBase<device_obj<LayerBase<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = LayerBase<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
        using device_obj_type = device_obj<Derived>;
    public:
        template<Scalar U>
        using MatrixND = DenseMatrix<U, MatrixOption::Col>;
        using ScalarType = Traits<device_obj<Derived>>::ScalarType;
        constexpr static bool IsTrain = ScalarType::isDiffable;
        constexpr static bool IsInfer = !IsTrain;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] auto forward(const auto& x) const { return Base::getDerived().forward(x); }
        auto reverse(const Derived& __restrict other) const noexcept;

        auto step(auto& optimizer) { return Base::getDerived().step(optimizer); }
        auto step() { return Base::getDerived().step(); }
        auto zero_grad() { return Base::getDerived().zero_grad(); }
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
