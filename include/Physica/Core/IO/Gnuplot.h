/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core {
    class Gnuplot {
    public:
        using VectorType = Vector<Scalar<Double, false>>;
        using DataArray = Utils::Array<Vector<Scalar<Double, false>>>;
    private:
        DataArray xDatas;
        DataArray yDatas;
    public:
        Gnuplot() = default;
        Gnuplot(const Gnuplot&) = default;
        Gnuplot(Gnuplot&&) noexcept = default;
        ~Gnuplot() = default;
        /* Operators */
        Gnuplot& operator=(const Gnuplot&) = default;
        Gnuplot& operator=(Gnuplot&&) noexcept = default;
        friend std::ostream& operator<<(std::ostream& os, const Gnuplot& poscar);
        friend std::istream& operator>>(std::istream& is, Gnuplot& poscar);
        /* Getters */
        [[nodiscard]] const DataArray& getXDatas() const noexcept { return xDatas; }
        [[nodiscard]] const DataArray& getYDatas() const noexcept { return yDatas; }
    };
}
