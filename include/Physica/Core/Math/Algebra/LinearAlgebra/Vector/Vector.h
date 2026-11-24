/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"

namespace Physica {
    template<class Derived, Scalar = Traits<Derived>::ScalarType> class RValueVector;

    template<class T>
    concept Vector = requires { typename std::remove_cvref_t<T>::ScalarType; } && []() consteval {
        using U = remove_codiff_t<std::remove_cvref_t<T>>;
        using S = U::ScalarType;
        return std::derived_from<U, RValueVector<U, S>>
            || std::derived_from<U, device_obj<RValueVector<typename remove_device_obj<U>::type, S>>>;
    }();
}
