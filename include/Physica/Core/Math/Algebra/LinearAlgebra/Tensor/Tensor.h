/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"

namespace Physica {
    template<class Derived, Scalar = Traits<Derived>::ScalarType> class RValueTensor;

    template<class T>
    concept Tensor = requires { typename std::remove_cvref_t<T>::ScalarType; } && []() consteval {
        using U = remove_codiff_t<std::remove_cvref_t<T>>;
        using S = U::ScalarType;
        return std::derived_from<U, RValueTensor<U, S>>
            || std::derived_from<U, device_obj<RValueTensor<typename remove_device_obj<U>::type, S>>>;
    }();
}
