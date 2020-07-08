/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_FUNCTIONDERIVATIVE_H
#define PHYSICA_FUNCTIONDERIVATIVE_H

#include <cstddef>

namespace Physica::Core {
    class Function;
    class FunctionTree;
    /*!
     * \class FunctionDerivative returns partial derivative of Function \f of variable \index.
     */
    class FunctionDerivative {
        const Function& f;
        size_t index;
    public:
        explicit FunctionDerivative(const Function& f, size_t index);

        [[nodiscard]] Function derivative() const;
    private:
        [[nodiscard]] FunctionTree derivativeTree(const FunctionTree& tree) const;
    };
}

#endif