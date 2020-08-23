/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_FUNCTIONDERIVATIVE_H
#define PHYSICA_FUNCTIONDERIVATIVE_H

#include <cstddef>
#include "TreeFunction.h"

namespace Physica::Core {
    /*!
     * \class FunctionDerivative returns partial derivative of TreeFunction \f of variable \index.
     */
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class FunctionDerivative {
        const TreeFunction<type, errorTrack>& f;
        size_t index;
    public:
        explicit FunctionDerivative(const TreeFunction<type, errorTrack>& f, size_t index);

        [[nodiscard]] TreeFunction<type, errorTrack> derivative() const;
    private:
        [[nodiscard]] TreeFunctionData<type, errorTrack> derivativeTree(const TreeFunctionData<type, errorTrack>& tree) const;
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/FunctionDerivativeImpl.h"

#endif