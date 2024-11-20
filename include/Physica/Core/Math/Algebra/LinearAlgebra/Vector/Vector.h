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

#include <concepts>
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    template<class Derived> class RValueVector;
    template<class Derived> class LValueVector;
    template<class Derived> class RSparseVector;

    template<class T>
    concept Vector = std::derived_from<T, RValueVector<T>> || std::derived_from<T, device_obj<RValueVector<T>>>;

    template<class T>
    concept LVector = std::derived_from<T, LValueVector<T>> || std::derived_from<T, device_obj<LValueVector<T>>>;

    template<class T>
    concept Sparse = std::derived_from<T, RSparseVector<T>>;
}
