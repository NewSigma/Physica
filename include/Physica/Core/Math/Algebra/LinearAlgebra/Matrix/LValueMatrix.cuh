/*
 * Copyright 2022 WeiBo He.
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

#include "RValueMatrix.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LValueMatrix<Derived>> : public device_obj<RValueMatrix<Derived>> {
        using Base = device_obj<RValueMatrix<Derived>>;
    public:
        using typename Base::ScalarType;
    public:
        /* Operators */
        [[nodiscard]] __device__ ScalarType& operator()(size_t row, size_t column) { return Base::getDerived()(row, column); }
        [[nodiscard]] __device__ const ScalarType& operator()(size_t row, size_t column) const { return Base::getDerived()(row, column); }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return (*this)(row, col); }
    };
}
