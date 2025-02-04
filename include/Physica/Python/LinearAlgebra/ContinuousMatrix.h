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

#include <pybind11/numpy.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/ContinuousMatrix.h"

namespace Physica::Core {
    template<class Derived>
    auto ContinuousMatrix<Derived>::toNumpy() const {
        static_assert(MatrixOption::isElementMatrix<Derived>(), "[Error]: Vector storage cannot cast to numpy");
        constexpr int Major = MatrixOption::getMajor<Derived>() == MatrixOption::Row
                            ? pybind11::array::c_style
                            : pybind11::array::f_style;
        using T = ScalarType::MachineType;
        return pybind11::array_t<T, Major>({Base::getRow(), Base::getCol()}, (const T*)data());
    }
}
