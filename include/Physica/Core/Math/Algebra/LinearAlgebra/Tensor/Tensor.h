/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"

namespace Physica {
    template<class Derived> class RValueTensor;

    namespace Internal {
        template<class T>
        concept TensorObj = std::derived_from<T, RValueTensor<T>>
                         || std::derived_from<T, device_obj<RValueTensor<typename remove_device_obj<T>::type>>>;
    }

    template<class T>
    concept Tensor = Internal::TensorObj<remove_codiff_t<std::remove_cvref_t<T>>>;
}
