/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_MATRIXCHAIN_H
#define PHYSICA_MATRIXCHAIN_H

#include <cstdlib>
#include <qglobal.h>
#include <memory>

namespace Physica::Core {
    class Matrix;
    class MatrixChain {
        Matrix** chain;
        size_t** price;
        size_t** point;
        size_t length;
    public:
        explicit MatrixChain(size_t length);
        ~MatrixChain();

        Matrix*& operator[](size_t i) { Q_ASSERT(i < length); return chain[i]; }
        std::unique_ptr<Matrix> solve();
    private:
        std::unique_ptr<Matrix> multiply(size_t from, size_t to);
    };
}

#endif
