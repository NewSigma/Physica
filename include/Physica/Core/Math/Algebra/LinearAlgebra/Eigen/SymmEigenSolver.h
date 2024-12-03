/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Decouplable.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.h"

namespace Physica::Core {
    /**
     * A = XTX^{-1}
     * where X is matrix of eigenvectors
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     */
    template<Scalar T, size_t Order = Dynamic>
    class SymmEigenSolver : public Decouplable {
        using This = SymmEigenSolver<T, Order>;
        using RealType = T::RealType;
        constexpr static bool isComplex = T::isComplex;
        static_assert(!isComplex, "[Error]: Complex matrix is not supported");
    public:
        using EigenvalueVector = DenseVector<RealType, Order>;
        using EigenvectorMatrix = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector, Order, Order>;
        using WorkingMatrix = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector, Order, Order>; // Optimize: Use tridiagonal matrix is better
    private:
        EigenvalueVector eigenvalues;
        EigenvectorMatrix eigenvectors;
        bool computeEigenvectors;
    public:
        SymmEigenSolver();
        SymmEigenSolver(size_t size);
        template<Matrix M>
        SymmEigenSolver(const M& source, bool computeEigenvectors_);
        SymmEigenSolver(const This&) = default;
        SymmEigenSolver(This&& solver) noexcept = default;
        ~SymmEigenSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        void compute(const M& source, bool computeEigenvectors_);
        void sort();
        template<class Compare>
        void sort(Compare comp);
        void resize(size_t size);
        void swap(This& __restrict solver) noexcept;
        /* Getters */
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] inline EigenvectorMatrix getEigenvectors() const noexcept;
    private:
        void stepQR(WorkingMatrix& working, size_t lower, size_t sub_order);
    };

    template<Scalar T, size_t Order>
    SymmEigenSolver<T, Order>::SymmEigenSolver() : eigenvalues(), eigenvectors(), computeEigenvectors(false) {}

    template<Scalar T, size_t Order>
    SymmEigenSolver<T, Order>::SymmEigenSolver(size_t size) : SymmEigenSolver() {
        resize(size);
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    SymmEigenSolver<T, Order>::SymmEigenSolver(const M& source, bool computeEigenvectors_)
            : SymmEigenSolver(source.getRow()) {
        compute(source, computeEigenvectors_);
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    void SymmEigenSolver<T, Order>::compute(const M& source, bool computeEigenvectors_) {
        assert(source.isSymm() && "[Error]: Bad symm matrix");
        assert(source.getRow() == eigenvalues.getLength() && "[Error]: Dimensions do not match");
        computeEigenvectors = computeEigenvectors_;
        [[unlikely]] if (source.getRow() == 1) {
            eigenvalues[0] = source.calc(0, 0);
            eigenvectors(0, 0) = T(1);
            return;
        }

        const RealType factor = abs_elem(source).max();
        if (factor < std::numeric_limits<T>::min()) {
            eigenvalues = RealType(0);
            return;
        }
        const RealType inv_factor = reciprocal(factor);
        const WorkingMatrix normalized = source * inv_factor; // Referenced from eigen, to avoid under/overflow in householder
        auto tridiagonal = Tridiagonalization<T, Order>(normalized);
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

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::sort() {
        sort([](T a, T b) noexcept { return a < b; });
    }

    template<Scalar T, size_t Order>
    template<class Compare>
    void SymmEigenSolver<T, Order>::sort(Compare comp) {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (comp(eigenvalues[j], eigenvalues[index_min]))
                    index_min = j;
            }

            const bool shouldSwap = i != index_min;
            if (!shouldSwap)
                continue;

            eigenvalues[i].swap(eigenvalues[index_min]);
            if (computeEigenvectors)
                eigenvectors.swap_col(i, index_min);
        }
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::resize(size_t size) {
        assert((Order == Dynamic || Order == size) && "[Error]: size is not consistent");
        eigenvalues.resize(size);
        eigenvectors.resize(size, size);
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::swap(This& __restrict solver) noexcept {
        assert(this != &solver && "[Error]: Self swap is likely a bug");
        eigenvalues.swap(solver.eigenvalues);
        std::swap(eigenvectors, solver.eigenvectors);
        std::swap(computeEigenvectors, solver.computeEigenvectors);
    }

    template<Scalar T, size_t Order>
    inline SymmEigenSolver<T, Order>::EigenvectorMatrix SymmEigenSolver<T, Order>::getEigenvectors() const noexcept {
        assert(computeEigenvectors && "[Error]: Eigenvectors are not ready");
        return eigenvectors;
    }
    /**
     * Use wilkinson shift
     */
    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::stepQR(WorkingMatrix& working, size_t lower, size_t sub_order) {
        Vector2D<T> buffer{};
        /* Init buffer */ {
            const auto subBlock = working.block(lower, sub_order, lower, sub_order);
            const RealType factor = (subBlock(sub_order - 2, sub_order - 2).real() - subBlock(sub_order - 1, sub_order - 1).real()) * T(0.5);
            const RealType factor2 = square(subBlock(sub_order - 1, sub_order - 2));
            const RealType factor3 = sqrt(square(factor) + factor2);
            const T shift = subBlock(sub_order - 1, sub_order - 1) - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3)); // TODO: why we introduce a divide operation
            buffer[0] = subBlock(0, 0) - shift;
            buffer[1] = subBlock(1, 0);
        }

        for (size_t i = 0; i < sub_order - 1; ++i) {
            const bool isShiftStep = i == 0;
            const size_t blockStart = lower + i + isShiftStep - 1;
            const size_t blockSize = ((i == sub_order - 2) ? 3 : 4) - isShiftStep;
            auto subBlock = working.block(blockStart, blockSize, blockStart, blockSize);
            if (!isShiftStep) {
                buffer[0] = subBlock(1, 0);
                buffer[1] = subBlock(2, 0);
            }

            const size_t index = !isShiftStep;
            auto givens_vec = givens(buffer, 0, 1);
            applyGivens(givens_vec, subBlock, index, index + 1);
            givens_vec[1] = -givens_vec[1];
            applyGivens(subBlock, givens_vec, index, index + 1);

            const T mean = (subBlock(1, index + index) + subBlock(index + 1, index)) * T(0.5);
            subBlock(index, index + 1) = subBlock(index + 1, index) = mean;

            if (computeEigenvectors)
                applyGivens(eigenvectors, givens_vec, lower + i, lower + i + 1);
            if (!isShiftStep)
                subBlock(2, 0) = subBlock(0, 2) = T(0);
        }
    }
}

namespace std {
    template<Physica::Core::Scalar T, size_t Order>
    inline void swap(Physica::Core::SymmEigenSolver<T, Order>& __restrict solver1,
                     Physica::Core::SymmEigenSolver<T, Order>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}
