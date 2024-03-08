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

#include <pybind11/pybind11.h>

namespace Physica::Python {
    class PhysicaObj {
        using This = PhysicaObj;
        py::bytes data;
    public:
        template<class T>
        PhysicaObj(T obj);
        PhysicaObj(const PhysicaObj&) = default;
        PhysicaObj(PhysicaObj&&) noexcept = default;
        ~PhysicaObj() = default;
        /* Operations */
        template<class T>
        [[nodiscard]] T& getDerived();
        template<class T>
        [[nodiscard]] const T& getDerived() const { return const_cast<This&>(*this).getDerived(); }
    };

    template<class T>
    PhysicaObj::PhysicaObj(T obj) {
        PlainStruct<T> buffer{};
        new (&buffer) T(std::move(obj));
        data = py::bytes((char*)(&buffer), sizeof(T));
    }

    template<class T>
    T& PhysicaObj::getDerived() {
        const std::string str(data);
        assert(sizeof(T) == str.size() && "[Error]: Type and bytes do not match");
        return reinterpret_cast<PlainStruct<T>*>(str.data())->getDerived();
    }
}
