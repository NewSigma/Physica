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
#include "Physica/Python/PythonMain.h"
#include "Physica/Python/Core/Algebra/Vector.h"
#include "Physica/Python/Core/MultiPrecision/Scalar.h"
#include "Physica/Python/Utils/Array.h"

PYBIND11_MODULE(PhysicaPython, m) {
    using namespace Physica::Python;
    m.doc() = "Python interface of Physica";
    Class<Scalar<Float>>::pybind(m);
    Class<Scalar<Double>>::pybind(m);
    Class<Physica::Utils::Array<Scalar<Double>>>::pybind(m);
    Class<Vector<Scalar<Double>>>::pybind(m);
}
