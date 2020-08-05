/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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