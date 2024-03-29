/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_ABSTRACTFUNCTION_H
#define PHYSICA_ABSTRACTFUNCTION_H

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    using namespace Utils;
    /*!
     * @class AbstractFunction contains common parts between @class TreeFunction and @class VectorFunction.
     */
    template<ScalarOption option, bool errorTrack>
    class AbstractFunction {
    protected:
        mutable Array<Scalar<option, errorTrack>> variables;
        mutable Array<Scalar<option, errorTrack>> constants;
    public:
        ~AbstractFunction() = default;
        /* Getters */
        [[nodiscard]] const Array<Scalar<option, errorTrack>>& getVariables() const { return variables; }
        [[nodiscard]] const Array<Scalar<option, errorTrack>>& getConstants() const { return constants; }
        [[nodiscard]] size_t getVariablePos(Scalar<option, errorTrack>* s) const;
        [[nodiscard]] size_t getConstantPos(Scalar<option, errorTrack>* s) const;
        /* Setters */
        inline void setConstant(Scalar<option, errorTrack> s, size_t index) const;
    protected:
        AbstractFunction(size_t variablesLength, size_t constantsLength);
        AbstractFunction(const AbstractFunction& f);
        AbstractFunction(AbstractFunction&& f) noexcept;
        /* Operators */
        AbstractFunction& operator=(const AbstractFunction& f);
        AbstractFunction& operator=(AbstractFunction&& f) noexcept;
        /* Setters */
        inline void setVariable(Scalar<option, errorTrack> s, size_t index) const;
    };
    /*!
     * Get the position of s in \variable, starts from 1.
     * return 0 if non-existent.
     */
    template<ScalarOption option, bool errorTrack>
    size_t AbstractFunction<option, errorTrack>::getVariablePos(Scalar<option, errorTrack>* s) const {
        const size_t distance = s - &variables[0];
        return distance < variables.getLength() ? distance + 1 : 0;
    }
    /*!
     * Get the position of s in \constant, starts from 1.
     * return 0 if non-existent.
     */
    template<ScalarOption option, bool errorTrack>
    size_t AbstractFunction<option, errorTrack>::getConstantPos(Scalar<option, errorTrack>* s) const {
        const size_t distance = s - &constants[0];
        return distance < variables.getLength() ? distance + 1 : 0;
    }

    template<ScalarOption option, bool errorTrack>
    inline void AbstractFunction<option, errorTrack>::setVariable(Scalar<option, errorTrack> s, size_t index) const {
        variables[index] = std::move(s);
    }

    template<ScalarOption option, bool errorTrack>
    inline void AbstractFunction<option, errorTrack>::setConstant(Scalar<option, errorTrack> s, size_t index) const {
        constants[index] = std::move(s);
    }
}

#include "AbstractFunctionImpl.h"

#endif
