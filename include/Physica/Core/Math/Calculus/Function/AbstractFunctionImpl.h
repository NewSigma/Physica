/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    AbstractFunction<option, errorTrack>::AbstractFunction(size_t variablesLength, size_t constantsLength)
            : variables(Array<Scalar<option, errorTrack>>(variablesLength))
            , constants(Array<Scalar<option, errorTrack>>(constantsLength)) {
        for(size_t i = 0; i < variablesLength; ++i)
            variables[i] = Scalar<option, errorTrack>(0);
        for(size_t i = 0; i < constantsLength; ++i)
            constants[i] = Scalar<option, errorTrack>(0);
    }

    template<ScalarOption option, bool errorTrack>
    AbstractFunction<option, errorTrack>::AbstractFunction(const AbstractFunction& f)
            : variables(f.variables), constants(f.constants) {}

    template<ScalarOption option, bool errorTrack>
    AbstractFunction<option, errorTrack>::AbstractFunction(AbstractFunction&& f) noexcept
            : variables(std::move(f.variables)), constants(std::move(f.constants)) {}
    /* Operators */
    template<ScalarOption option, bool errorTrack>
    AbstractFunction<option, errorTrack>& AbstractFunction<option, errorTrack>::operator=(const AbstractFunction& f) {
        if(this != &f) {
            variables = f.variables;
            constants = f.constants;
        }
        return *this;
    }
    template<ScalarOption option, bool errorTrack>
    AbstractFunction<option, errorTrack>& AbstractFunction<option, errorTrack>::operator=(AbstractFunction&& f) noexcept {
        variables = std::move(f.variables);
        constants = std::move(f.constants);
        return *this;
    }
}
