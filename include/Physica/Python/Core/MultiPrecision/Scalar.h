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

#include "Physica/Core/MultiPrecision/Scalar.h"

using namespace Physica::Core;

namespace Physica::Python {
    template<ScalarOption option>
    class Class<Scalar<option>> {
        using BindType = Scalar<option>;
        constexpr static const char* BindName = option == Float ? "Scalar<Float>" : "Scalar<Double>";
    public:
        inline static void pybind(::pybind11::module_& m) noexcept {
            py::class_<BindType>(m, BindName)
                .def("__float__", &BindType::getTrivial);
        }
    };
}
