/*
 * Copyright 2022 WeiBo He.
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
#pragma once

#include "Matrix/MatrixDecomposition/Tridiagonalization.h"
#include "Matrix/MatrixDecomposition/Schur.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     */
    template<class MatrixType>
    class SymmEigenSolver : public AbstractSchur {
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;
        constexpr static bool isComplex = ScalarType::isComplex;
        static_assert(!isComplex, "Complex matrix is not supported");
    public:
        using EigenvalueVector = Vector<RealType, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
        using EigenvectorMatrix = DenseMatrix<ComplexScalar<RealType>,
                                              DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                              MatrixType::RowAtCompile,
                                              MatrixType::RowAtCompile>;
        using WorkingMatrix = DenseMatrix<ScalarType,
                                          DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                          MatrixType::RowAtCompile,
                                          MatrixType::RowAtCompile>; //Optimize: Use tridiagonal matrix is better
    private:
        EigenvalueVector eigenvalues;
        EigenvectorMatrix eigenvectors;
        bool computeEigenvectors;
    public:
        SymmEigenSolver();
        SymmEigenSolver(size_t size);
        template<class OtherMatrix>
        SymmEigenSolver(const RValueMatrix<OtherMatrix>& source, bool computeEigenvectors_);
        SymmEigenSolver(const SymmEigenSolver&) = default;
        SymmEigenSolver(SymmEigenSolver&& solver) noexcept = default;
        ~SymmEigenSolver() = default;
        /* Operators */
        SymmEigenSolver& operator=(SymmEigenSolver solver) noexcept;
        /* Operations */
        template<class OtherMatrix>
        void compute(const RValueMatrix<OtherMatrix>& source, bool computeEigenvectors_);
        void sort();
        /* Getters */
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] EigenvectorMatrix getEigenvectors() const noexcept { return eigenvectors; }
        /* Helpers */
        void swap(SymmEigenSolver& solver) noexcept;
    private:
        static size_t activeWindowLower(WorkingMatrix& mat, size_t upper);
        void stepQR(WorkingMatrix& working, size_t lower, size_t sub_order);
    };

    template<class MatrixType>
    SymmEigenSolver<MatrixType>::SymmEigenSolver() : eigenvalues(), eigenvectors(), computeEigenvectors(false) {}

    template<class MatrixType>
    SymmEigenSolver<MatrixType>::SymmEigenSolver(size_t size)
            : eigenvalues(size)
            , eigenvectors(size, size) {}

    template<class MatrixType>
    template<class OtherMatrix>
    SymmEigenSolver<MatrixType>::SymmEigenSolver(const RValueMatrix<OtherMatrix>& source, bool computeEigenvectors_)
            : eigenvalues(source.getRow())
            , eigenvectors(source.getRow(), source.getRow())
            , computeEigenvectors(computeEigenvectors_) {
        compute(source, computeEigenvectors);
    }

    template<class MatrixType>
    SymmEigenSolver<MatrixType>& SymmEigenSolver<MatrixType>::operator=(SymmEigenSolver<MatrixType> solver) noexcept {
        swap(solver);
        return *this;
    }

    template<class MatrixType>
    template<class OtherMatrix>
    void SymmEigenSolver<MatrixType>::compute(const RValueMatrix<OtherMatrix>& source, bool computeEigenvectors_) {
        assert(source.getRow() == source.getColumn());
        assert(source.getRow() == eigenvalues.getLength());
        computeEigenvectors = computeEigenvectors_;
        
        typename MatrixType::RealMatrix buffer = abs(source);
        const RealType factor = buffer.max();
        if (factor < std::numeric_limits<ScalarType>::min()) {
            eigenvalues = RealType::Zero();
            return;
        }
        const RealType inv_factor = reciprocal(factor);
        const MatrixType normalized = source * inv_factor; //Referenced from eigen, to avoid under/overflow in householder
        auto tridiagonal = Tridiagonalization<MatrixType>(normalized);
        WorkingMatrix working = tridiagonal.getMatrixT();
        if (computeEigenvectors)
            eigenvectors = tridiagonal.getMatrixQ();

        const size_t order = working.getRow();
        size_t upper = order - 1;
        size_t iter = 0;
        size_t total_iter = 0;
        const size_t max_iter = AbstractSchur::maxItePerCol * order;
        while (1 <= upper && upper < order) {
            const size_t lower = AbstractSchur::activeWindowLower(working, upper);
            if (lower == upper) {
                upper -= 1;
                iter = 0;
            }
            else {
                const size_t sub_order = upper - lower + 1;
                stepQR(working, lower, sub_order);
                ++iter;
                ++total_iter;
            }

            if (total_iter == max_iter)
                throw BadConvergenceException();
        }

        for (size_t i = 0; i < order; ++i)
            eigenvalues[i] = working(i, i) * factor;
    }

    template<class MatrixType>
    void SymmEigenSolver<MatrixType>::sort() {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (eigenvalues[j] < eigenvalues[index_min])
                    index_min = j;
            }
            eigenvalues[i].swap(eigenvalues[index_min]);
            if (computeEigenvectors)
                eigenvectors[i].swap(eigenvectors[index_min]);
        }
    }

    template<class MatrixType>
    void SymmEigenSolver<MatrixType>::swap(SymmEigenSolver<MatrixType>& solver) noexcept {
        eigenvalues.swap(solver.eigenvalues);
        std::swap(eigenvectors, solver.eigenvectors);
        std::swap(computeEigenvectors, solver.computeEigenvectors);
    }

    template<class MatrixType>
    void SymmEigenSolver<MatrixType>::stepQR(WorkingMatrix& working, size_t lower, size_t sub_order) {
        auto subBlock = working.block(lower, sub_order, lower, sub_order);
        const RealType factor = (subBlock(sub_order - 2, sub_order - 2).getReal() - subBlock(sub_order - 1, sub_order - 1).getReal()) * 0.5;
        const RealType factor2 = square(subBlock(sub_order - 1, sub_order - 2));
        const RealType factor3 = sqrt(square(factor) + factor2);
        const ScalarType shift = subBlock(sub_order - 1, sub_order - 1) - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3));
        
        Vector<ScalarType, 2> buffer{subBlock(0, 0) - shift, subBlock(1, 0)};
        auto givens_vec = givens(buffer, 0, 1);
        applyGivens(givens_vec, subBlock, 0, 1);
        givens_vec[1].toOpposite();
        applyGivens(subBlock, givens_vec, 0, 1);
        if (computeEigenvectors)
            applyGivens(eigenvectors, givens_vec, 0, 1);
        for (size_t i = 1; i < sub_order - 1; ++i) {
            buffer[0] = subBlock(i, i - 1);
            buffer[1] = subBlock(i + 1, i - 1);
            givens_vec = givens(buffer, 0, 1);
            applyGivens(givens_vec, subBlock, i, i + 1);
            givens_vec[1].toOpposite();
            applyGivens(subBlock, givens_vec, i, i + 1);
            if (computeEigenvectors)
                applyGivens(eigenvectors, givens_vec, i, i + 1);
        }
    }

    template<class MatrixType>
    inline void swap(Physica::Core::SymmEigenSolver<MatrixType>& solver1,
                     Physica::Core::SymmEigenSolver<MatrixType>& solver2) noexcept {
        solver1.swap(solver2);
    }
}
