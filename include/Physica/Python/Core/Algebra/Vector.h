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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

using namespace Physica::Core;

namespace Physica::Python {
    template<class T, size_t Length, size_t MaxLength, class Allocator>
    class Class<Vector<T, Length, MaxLength, Allocator>> {
        using BindType = Vector<T, Length, MaxLength, Allocator>;
        using Storage = typename BindType::Storage;
        using TrivialType = typename T::TrivialType;
        constexpr static bool isFloat = T::Option == Float;
        constexpr static const char* BindName = isFloat ? "Vector<Float>" : "Vector<Double>";
    public:
        inline static void pybind(::pybind11::module_& m) noexcept {
            py::class_<BindType, Storage>(m, BindName, py::buffer_protocol())
                .def_buffer([](BindType& v) -> py::buffer_info {
                    return py::buffer_info(
                        v.data(),
                        sizeof(double),
                        py::format_descriptor<double>::format(),
                        1,
                        {v.getLength()},
                        {sizeof(double)}
                    );
                })
                .def("calc", &BindType::calc);
        }
    };
}

