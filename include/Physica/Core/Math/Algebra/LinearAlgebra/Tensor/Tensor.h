/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica {
    template<class Derived> class RValueTensor;

    namespace Internal {
        template<class T>
        concept TensorObj = std::derived_from<T, RValueTensor<T>>
                         || std::derived_from<T, device_obj<RValueTensor<typename remove_device_obj<T>::Type>>>;
    }

    template<class T>
    concept Tensor = Internal::TensorObj<std::remove_cvref_t<T>> || Internal::TensorObj<typename remove_codiff<std::remove_cvref_t<T>>::Type>;

    using Index2D = Array<size_t, 2>;
    using Index3D = Array<size_t, 3>;
    using Index4D = Array<size_t, 4>;
    using Index5D = Array<size_t, 5>;
    using IndexND = Array<size_t>;
}
