/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LUDECOMPOSITION_H
#define PHYSICA_LUDECOMPOSITION_H

#include <cstdlib>
#include "Scalar.h"
#include "Matrix.h"

namespace Physica::Core {
    class SquareMatrix;
    class LUDecomposition {
        Matrix* matrix;
        size_t* biasOrder;
    public:
        explicit LUDecomposition(SquareMatrix& square);
        explicit LUDecomposition(const SquareMatrix& square);
        ~LUDecomposition();
        /* Getters */
        [[nodiscard]] Matrix* release() { auto temp = matrix; matrix = nullptr; return temp; }
        [[nodiscard]] const Matrix& getMatrix() { return *matrix; }
        [[nodiscard]] const size_t* getOrder() { return biasOrder; }
    private:
        static void decompositionColumn(Matrix& m, size_t column);
    };
}

#endif
