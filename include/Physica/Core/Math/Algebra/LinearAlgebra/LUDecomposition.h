/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LUDECOMPOSITION_H
#define PHYSICA_LUDECOMPOSITION_H

#include <cstdlib>
#include "Matrix/Matrix.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    class LUDecomposition {
        Matrix<T, type, maxRow, maxColumn> matrix;
        size_t* biasOrder;
    public:
        explicit LUDecomposition(const Matrix<T, type, maxRow, maxColumn>& m);
        explicit LUDecomposition(Matrix<T, type, maxRow, maxColumn>&& m) noexcept;
        LUDecomposition(const LUDecomposition& l);
        LUDecomposition(LUDecomposition&& l) noexcept;
        ~LUDecomposition();
        /* Operators */
        LUDecomposition& operator=(const LUDecomposition& l);
        LUDecomposition& operator=(LUDecomposition&& l) noexcept;
        /* Getters */
        [[nodiscard]] Matrix<T, type, maxRow, maxColumn>&& release() { return std::move(matrix); }
        [[nodiscard]] const Matrix<T, type, maxRow, maxColumn>& getMatrix() { return matrix; }
        [[nodiscard]] const size_t* getOrder() { return biasOrder; }
    private:
        void decompositionColumn(size_t column);
    };
}

#include "LUDecompositionImpl.h"

#endif
