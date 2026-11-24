/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "LayerMixin.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"

namespace Physica {
    template<class Derived>
    class device_obj<LayerMixin<Derived>> : public CRTP<device_obj<LayerMixin<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = LayerMixin<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTP<This>;
        using device_obj_type = device_obj<Derived>;
    public:
        template<Scalar U>
        using MatrixND = DenseMatrix<U, MatrixMajor::Col>;
        using ScalarType = Traits<device_obj<Derived>>::ScalarType;
    protected:
        using T = ScalarType;
        using Tv = T::ValueType;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] auto forward(const auto& x) const { return Base::getDerived().forward(x); }
        auto reverse(const Derived& __restrict other) const noexcept;

        auto step(auto& optimizer) { return Base::getDerived().step(optimizer); }
        auto step() { return Base::getDerived().step(); }
        auto zero_grad() { return Base::getDerived().zero_grad(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ consteval static bool isTraining() noexcept { return T::isDiffable(); }
        [[nodiscard]] __host__ __device__ consteval static bool isInfering() noexcept { return !isTraining(); }
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
    class Traits<device_obj<LayerMixin<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}
