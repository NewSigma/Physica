/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTFUNCTIONIMPL_H
#define PHYSICA_ABSTRACTFUNCTIONIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    AbstractFunction<type, errorTrack>::AbstractFunction(size_t variablesLength, size_t constantsLength)
            : variables(CStyleArray<Scalar<type, errorTrack>, Dynamic>(variablesLength))
            , constants(CStyleArray<Scalar<type, errorTrack>, Dynamic>(constantsLength)) {}

    template<ScalarType type, bool errorTrack>
    AbstractFunction<type, errorTrack>::AbstractFunction(const AbstractFunction& f)
            : variables(f.variables), constants(f.constants) {}
    template<ScalarType type, bool errorTrack>
    AbstractFunction<type, errorTrack>::AbstractFunction(AbstractFunction&& f) noexcept
            : variables(std::move(f.variables)), constants(std::move(f.constants)) {}
    /* Operators */
    template<ScalarType type, bool errorTrack>
    AbstractFunction<type, errorTrack>& AbstractFunction<type, errorTrack>::operator=(const AbstractFunction& f) {
        if(this != &f) {
            variables = f.variables;
            constants = f.constants;
        }
        return *this;
    }
    template<ScalarType type, bool errorTrack>
    AbstractFunction<type, errorTrack>& AbstractFunction<type, errorTrack>::operator=(AbstractFunction&& f) noexcept {
        variables = std::move(f.variables);
        constants = std::move(f.constants);
        return *this;
    }
}

#endif
