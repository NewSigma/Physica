/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_MATRIXCHAIN_H
#define PHYSICA_MATRIXCHAIN_H

#include <cstdlib>
#include <qglobal.h>
#include <memory>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    class MatrixChain {
        Matrix<T, type, maxRow, maxColumn>** chain;
        size_t** price;
        size_t** point;
        size_t length;
    public:
        explicit MatrixChain(size_t length);
        ~MatrixChain();

        Matrix<T, type, maxRow, maxColumn>*& operator[](size_t i) { Q_ASSERT(i < length); return chain[i]; }
        Matrix<T, type, Dynamic, Dynamic> solve();
    private:
        Matrix<T, type, Dynamic, Dynamic> multiply(size_t from, size_t to);
    };
}

#endif
