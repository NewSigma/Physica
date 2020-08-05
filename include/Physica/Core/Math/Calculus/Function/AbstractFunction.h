/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTFUNCTION_H
#define PHYSICA_ABSTRACTFUNCTION_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/Container/CStyleArray/CStyleArray.h"

namespace Physica::Core {
    /*!
     * @class AbstractFunction contains common parts between @class TreeFunction and @class VectorFunction.
     */
    template<ScalarType type, bool errorTrack>
    class AbstractFunction {
    protected:
        CStyleArray<Scalar<type, errorTrack>, Dynamic> variables;
        CStyleArray<Scalar<type, errorTrack>, Dynamic> constants;
    public:
        ~AbstractFunction() = default;
    protected:
        AbstractFunction(size_t variablesLength, size_t constantsLength);
        AbstractFunction(const AbstractFunction& f);
        AbstractFunction(AbstractFunction&& f) noexcept;
        /* Operators */
        AbstractFunction& operator=(const AbstractFunction& f);
        AbstractFunction& operator=(AbstractFunction&& f) noexcept;
        /* Getters */
        [[nodiscard]] const CStyleArray<Scalar<type, errorTrack>, Dynamic>& getVariables() { return variables; }
        [[nodiscard]] const CStyleArray<Scalar<type, errorTrack>, Dynamic>& getConstants() { return constants; }
        [[nodiscard]] size_t getVariablePos(Scalar<type, errorTrack>* s) const;
        [[nodiscard]] size_t getConstantPos(Scalar<type, errorTrack>* s) const;
    };
    /*!
     * Get the position of s in \variable, starts from 1.
     * return 0 if non-existent.
     */
    template<ScalarType type, bool errorTrack>
    size_t AbstractFunction<type, errorTrack>::getVariablePos(Scalar<type, errorTrack>* s) const {
        const size_t distance = s - &variables[0];
        return distance < variables.getLength() ? distance + 1 : 0;
    }
    /*!
     * Get the position of s in \constant, starts from 1.
     * return 0 if non-existent.
     */
    template<ScalarType type, bool errorTrack>
    size_t AbstractFunction<type, errorTrack>::getConstantPos(Scalar<type, errorTrack>* s) const {
        const size_t distance = s - &constants[0];
        return distance < variables.getLength() ? distance + 1 : 0;
    }
}

#include "AbstractFunctionImpl.h"

#endif
