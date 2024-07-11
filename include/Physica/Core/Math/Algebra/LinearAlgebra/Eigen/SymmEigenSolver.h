/*
 * Copyright 2022-2024 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Schur.h"

namespace Physica::Core {
    /**
     * A = XTX^{-1}
     * where X is matrix of eigenvectors
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     */
    template<class ScalarType, size_t Order = Dynamic>
    class SymmEigenSolver : public Decouplable {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        using This = SymmEigenSolver<ScalarType, Order>;
        using RealType = typename ScalarType::RealType;
        constexpr static bool isComplex = ScalarType::isComplex;
        static_assert(!isComplex, "[Error]: Complex matrix is not supported");
    public:
        using EigenvalueVector = Vector<RealType, Order>;
        using EigenvectorMatrix = DenseMatrix<ComplexScalar<RealType>, MatrixOption::Column | MatrixOption::Vector, Order, Order>;
        using WorkingMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector, Order, Order>; //Optimize: Use tridiagonal matrix is better
    private:
        EigenvalueVector eigenvalues;
        EigenvectorMatrix eigenvectors;
        bool computeEigenvectors;
    public:
        SymmEigenSolver();
        SymmEigenSolver(size_t size);
        template<class MatrixType>
        SymmEigenSolver(const RValueMatrix<MatrixType>& source, bool computeEigenvectors_);
        SymmEigenSolver(const This&) = default;
        SymmEigenSolver(This&& solver) noexcept = default;
        ~SymmEigenSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class MatrixType>
        void compute(const RValueMatrix<MatrixType>& source, bool computeEigenvectors_);
        void sort();
        void resize(size_t size);
        void swap(This& __restrict solver) noexcept;
        /* Getters */
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] inline EigenvectorMatrix getEigenvectors() const noexcept;
    private:
        void stepQR(WorkingMatrix& working, size_t lower, size_t sub_order);
    };

    template<class ScalarType, size_t Order>
    SymmEigenSolver<ScalarType, Order>::SymmEigenSolver() : eigenvalues(), eigenvectors(), computeEigenvectors(false) {}

    template<class ScalarType, size_t Order>
    SymmEigenSolver<ScalarType, Order>::SymmEigenSolver(size_t size) : SymmEigenSolver() {
        resize(size);
    }

    template<class ScalarType, size_t Order>
    template<class MatrixType>
    SymmEigenSolver<ScalarType, Order>::SymmEigenSolver(const RValueMatrix<MatrixType>& source, bool computeEigenvectors_)
            : SymmEigenSolver(source.getRow()) {
        compute(source, computeEigenvectors_);
    }

    template<class ScalarType, size_t Order>
    template<class MatrixType>
    void SymmEigenSolver<ScalarType, Order>::compute(const RValueMatrix<MatrixType>& source, bool computeEigenvectors_) {
        assert(source.getRow() == source.getColumn() && "[Error]: Square matrix is required");
        assert(source.getRow() == eigenvalues.getLength());
        computeEigenvectors = computeEigenvectors_;
        [[unlikely]] if (source.getRow() == 1) {
            eigenvalues[0] = source.calc(0, 0);
            eigenvectors(0, 0) = ScalarType(1);
            return;
        }

        const RealType factor = abs_elem(source).max();
        if (factor < std::numeric_limits<ScalarType>::min()) {
            eigenvalues = RealType(0);
            return;
        }
        const RealType inv_factor = reciprocal(factor);
        const WorkingMatrix normalized = source * inv_factor; //Referenced from eigen, to avoid under/overflow in householder
        auto tridiagonal = Tridiagonalization<ScalarType, Order>(normalized);
        WorkingMatrix working = tridiagonal.getMatrixT();
        if (computeEigenvectors)
            eigenvectors = tridiagonal.getMatrixQ();

        const size_t order = working.getRow();
        size_t upper = order - 1;
        size_t total_iter = 0;
        const size_t max_iter = Decouplable::maxItePerCol * order;
        while (1 <= upper && upper < order) {
            const size_t lower = Decouplable::activeWindowDownDiag(working, upper);
            if (lower == upper) {
                upper -= 1;
            }
            else {
                const size_t sub_order = upper - lower + 1;
                stepQR(working, lower, sub_order);
                ++total_iter;
            }

            if (total_iter == max_iter)
                throw BadConvergenceException("Exceed max iteration of SymmEigenSolver");
        }

        for (size_t i = 0; i < order; ++i)
            eigenvalues[i] = working(i, i) * factor;
    }

    template<class ScalarType, size_t Order>
    void SymmEigenSolver<ScalarType, Order>::sort() {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (eigenvalues[j] < eigenvalues[index_min])
                    index_min = j;
            }

            const bool shouldSwap = i != index_min;
            if (!shouldSwap)
                continue;

            eigenvalues[i].swap(eigenvalues[index_min]);
            if (computeEigenvectors)
                eigenvectors.asArray()[i].swap(eigenvectors.asArray()[index_min]);
        }
    }

    template<class ScalarType, size_t Order>
    void SymmEigenSolver<ScalarType, Order>::resize(size_t size) {
        assert((Order == Dynamic || Order == size) && "[Error]: size is not consistent");
        eigenvalues.resize(size);
        eigenvectors.resize(size, size);
    }

    template<class ScalarType, size_t Order>
    void SymmEigenSolver<ScalarType, Order>::swap(This& __restrict solver) noexcept {
        assert(this != &solver && "[Error]: Self swap is likely a bug");
        eigenvalues.swap(solver.eigenvalues);
        std::swap(eigenvectors, solver.eigenvectors);
        std::swap(computeEigenvectors, solver.computeEigenvectors);
    }

    template<class ScalarType, size_t Order>
    inline typename SymmEigenSolver<ScalarType, Order>::EigenvectorMatrix SymmEigenSolver<ScalarType, Order>::getEigenvectors() const noexcept {
        assert(computeEigenvectors && "[Error]: Eigenvectors are not ready");
        return eigenvectors;
    }
    /**
     * Use wilkinson shift
     */
    template<class ScalarType, size_t Order>
    void SymmEigenSolver<ScalarType, Order>::stepQR(WorkingMatrix& working, size_t lower, size_t sub_order) {
        auto subBlock = working.block(lower, sub_order, lower, sub_order);
        const RealType factor = (subBlock(sub_order - 2, sub_order - 2).getReal() - subBlock(sub_order - 1, sub_order - 1).getReal()) * 0.5;
        const RealType factor2 = square(subBlock(sub_order - 1, sub_order - 2));
        const RealType factor3 = sqrt(square(factor) + factor2);
        const ScalarType shift = subBlock(sub_order - 1, sub_order - 1) - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3)); //TODO: why we introduce a divide operation
        
        Vector<ScalarType, 2> buffer{subBlock(0, 0) - shift, subBlock(1, 0)};
        auto givens_vec = givens(buffer, 0, 1);
        applyGivens(givens_vec, subBlock, 0, 1);
        givens_vec[1].toOpposite();
        applyGivens(subBlock, givens_vec, 0, 1);
        if (computeEigenvectors)
            applyGivens(eigenvectors, givens_vec, lower, lower + 1);
        for (size_t i = 1; i < sub_order - 1; ++i) {
            buffer[0] = subBlock(i, i - 1);
            buffer[1] = subBlock(i + 1, i - 1);
            givens_vec = givens(buffer, 0, 1);
            applyGivens(givens_vec, subBlock, i, i + 1);
            givens_vec[1].toOpposite();
            applyGivens(subBlock, givens_vec, i, i + 1);
            if (computeEigenvectors)
                applyGivens(eigenvectors, givens_vec, lower + i, lower + i + 1);
        }
    }
}

namespace std {
    template<class ScalarType, size_t Order>
    inline void swap(Physica::Core::SymmEigenSolver<ScalarType, Order>& __restrict solver1,
                     Physica::Core::SymmEigenSolver<ScalarType, Order>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}
