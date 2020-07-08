/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LINEAREQUATIONS_H
#define PHYSICA_LINEAREQUATIONS_H

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<class T = MultiScalar, MatrixType type = Column, size_t maxRow = Dynamic, size_t maxColumn = Dynamic>
    class LinearEquations {
    public:
        enum LinearEquationsMethod {
            GaussJordanPartial,
            GaussJordanComplete,
            GaussEliminationPartial,
            GaussEliminationComplete,
            LUMethod
        };
    private:
        Matrix<T, type, maxRow, maxColumn> matrix;
    public:
        explicit LinearEquations(const Matrix<T, type, maxRow, maxColumn>& m);
        explicit LinearEquations(Matrix<T, type, maxRow, maxColumn>&& m) noexcept;
        LinearEquations(const LinearEquations& l) = default;
        LinearEquations(LinearEquations&& l) noexcept;
        ~LinearEquations() = default;
        /* Operators */
        LinearEquations& operator=(const LinearEquations& l) = default;
        LinearEquations& operator=(LinearEquations&& l) noexcept;
        /* Operations */
        const typename Matrix<T, type, maxRow, maxColumn>::VectorType& solve(LinearEquationsMethod method);
        /* Helpers */
        [[nodiscard]] Matrix<T, type, maxRow, maxColumn>&& release() noexcept { return std::move(matrix); }
        /* Getters */
        [[nodiscard]] const Matrix<T, type, maxRow, maxColumn>& getMatrix() const noexcept { return matrix; }
    };
}

#include "LinearEquationsImpl.h"

#endif