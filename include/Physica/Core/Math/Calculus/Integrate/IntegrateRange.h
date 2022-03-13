/*
 * Copyright 2020-2022 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, size_t dim>
    class IntegrateRange {
        //We consider dimension which is larger than 3 is unsuitable to allocate domain data on stack.
        static_assert(dim <= 3, "Dimension larger than 3 must be set Dynamic.");
    public:
        using VectorType = Vector<ScalarType, dim>;
    private:
        VectorType range_from;
        VectorType range_to;
    public:
        IntegrateRange(VectorType range_from_, VectorType range_to_);
        /* Getters */
        [[nodiscard]] const VectorType& from() const { return range_from; }
        [[nodiscard]] const VectorType& to() const { return range_to; }
    };

    template<class ScalarType, size_t dim>
    IntegrateRange<ScalarType, dim>::IntegrateRange(VectorType range_from_, VectorType range_to_)
            : range_from(std::move(range_from_)), range_to(std::move(range_to_)) {}
}
