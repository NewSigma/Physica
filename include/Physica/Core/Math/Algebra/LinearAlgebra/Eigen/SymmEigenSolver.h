/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Decouplable.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Tridiagonalization.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
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
        constexpr static bool isComplex = T::isComplex;
        static_assert(!isComplex, "[Error]: Complex matrix is not supported");

        using Tr = T::RealType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        using EigenvalueVector = DenseVector<Tr, Order>;
        using EigenvectorMatrix = DenseMatrix<T, MatrixOption::Col, Order, Order>;
        using WorkingMatrix = DenseMatrix<T, MatrixOption::Col, Order, Order>; // Optimize: Use tridiagonal matrix is better
    private:
        EigenvalueVector eigenvalues;
        EigenvectorMatrix eigenvectors;
    public:
        SymmEigenSolver() = default;
        SymmEigenSolver(size_t size, bool needEigenvectors);
        SymmEigenSolver(const Matrix auto& source, bool needEigenvectors);
        SymmEigenSolver(const This&) = default;
        SymmEigenSolver(This&& solver) noexcept = default;
        ~SymmEigenSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void compute_mkl(const Matrix auto& source);
        void compute_base(const Matrix auto& source);

        void sort();
        void sort(std::invocable<Tr, Tr> auto comp);
        void resize(size_t size);
        void resize(size_t size, bool needEigenvectors);
        void swap(This& __restrict solver) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getEigenvalues(this auto&& self) noexcept;
        [[nodiscard]] const EigenvectorMatrix& getEigenvectors() const noexcept;
        [[nodiscard]] bool getNeedEigenvectors() const noexcept;
    private:
        void stepQR(WorkingMatrix& working, size_t lower, size_t sub_order);
        /* Static members */
        static bool defaultComp(T a, T b) noexcept;
    };

    template<Scalar T, size_t Order>
    SymmEigenSolver<T, Order>::SymmEigenSolver(size_t size, bool needEigenvectors) {
        resize(size, needEigenvectors);
    }

    template<Scalar T, size_t Order>
    SymmEigenSolver<T, Order>::SymmEigenSolver(const Matrix auto& source, bool needEigenvectors)
            : SymmEigenSolver(source.getRow(), needEigenvectors) {
        compute(source);
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::compute(const Matrix auto& source) {
        if constexpr (HasMKL() && !Diffable<T>)
            compute_mkl(source);
        else
            compute_base(source);
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::compute_base(const Matrix auto& source) {
        assert(source.isSymm() && "[Error]: Bad symm matrix");
        assert(source.getRow() == eigenvalues.getLength() && "[Error]: Dimensions do not match");
        if (source.getRow() == 1) [[unlikely]] {
            eigenvalues[0] = source.calc(0, 0);
            if (getNeedEigenvectors())
                eigenvectors[0, 0] = T(1);
            return;
        }

        const Tr factor = abs_elem(source).max();
        if (factor.isSubNormal()) {
            eigenvalues = Tr(0);
            return;
        }
        const Tr inv_factor = reciprocal(factor);
        const WorkingMatrix normalized = source * inv_factor; // Referenced from eigen, to avoid under/overflow in householder
        auto tridiagonal = Tridiagonalization<T, Order>(normalized);
        WorkingMatrix working = tridiagonal.getMatrixT();
        if (getNeedEigenvectors())
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
            eigenvalues[i] = working[i, i] * factor;
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::sort() {
        sort(defaultComp);
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::sort(std::invocable<Tr, Tr> auto comp) {
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
            if (getNeedEigenvectors())
                eigenvectors.swap_col(i, index_min);
        }
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::resize(size_t size) {
        resize(size, getNeedEigenvectors());
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::resize(size_t size, bool needEigenvectors) {
        assert((Order == Dynamic || Order == size) && "[Error]: size is not consistent");
        eigenvalues.resize(size);
        if (needEigenvectors)
            eigenvectors.resize(size, size);
    }

    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::swap(This& __restrict solver) noexcept {
        assert(this != &solver && "[Error]: Self swap is likely a bug");
        eigenvalues.swap(solver.eigenvalues);
        eigenvectors.swap(solver.eigenvectors);
    }

    template<Scalar T, size_t Order>
    auto&& SymmEigenSolver<T, Order>::getEigenvalues(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).eigenvalues;
    }

    template<Scalar T, size_t Order>
    auto SymmEigenSolver<T, Order>::getEigenvectors() const noexcept -> const EigenvectorMatrix& {
        assert(getNeedEigenvectors() && "[Error]: Eigenvectors are not ready");
        return eigenvectors;
    }

    template<Scalar T, size_t Order>
    bool SymmEigenSolver<T, Order>::getNeedEigenvectors() const noexcept {
        return !eigenvectors.empty();
    }
    /**
     * Use wilkinson shift
     */
    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::stepQR(WorkingMatrix& working, size_t lower, size_t sub_order) {
        Vector2D<T> buffer{};
        /* Init buffer */ {
            const auto subBlock = working.block(lower, sub_order, lower, sub_order);
            const Tr factor = (subBlock[sub_order - 2, sub_order - 2].real() - subBlock[sub_order - 1, sub_order - 1].real()) * T(0.5);
            const Tr factor2 = square(subBlock[sub_order - 1, sub_order - 2]);
            const Tr factor3 = sqrt(fma(factor, factor, factor2));
            const T shift = subBlock[sub_order - 1, sub_order - 1] - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3)); // TODO: why we introduce a divide operation
            buffer[0] = subBlock[0, 0] - shift;
            buffer[1] = subBlock[1, 0];
        }

        for (size_t i = 0; i < sub_order - 1; ++i) {
            const bool isShiftStep = i == 0;
            const size_t blockStart = lower + i + isShiftStep - 1;
            const size_t blockSize = ((i == sub_order - 2) ? 3 : 4) - isShiftStep;
            auto subBlock = working.block(blockStart, blockSize, blockStart, blockSize);
            if (!isShiftStep) {
                buffer[0] = subBlock[1, 0];
                buffer[1] = subBlock[2, 0];
            }

            const size_t index = !isShiftStep;
            auto givens_vec = givens(buffer, 0, 1);
            applyGivens(givens_vec, subBlock, index, index + 1);
            givens_vec[1] = -givens_vec[1];
            applyGivens(subBlock, givens_vec, index, index + 1);

            const T mean = (subBlock[index, index + 1] + subBlock[index + 1, index]) * T(0.5);
            subBlock[index, index + 1] = subBlock[index + 1, index] = mean;

            if (getNeedEigenvectors())
                applyGivens(eigenvectors, givens_vec, lower + i, lower + i + 1);
            if (!isShiftStep)
                subBlock[2, 0] = subBlock[0, 2] = T(0);
        }
    }

    template<Scalar T, size_t Order>
    bool SymmEigenSolver<T, Order>::defaultComp(T a, T b) noexcept {
        return a < b;
    }

    template<Scalar T, size_t Order>
    void swap(SymmEigenSolver<T, Order>& __restrict solver1, SymmEigenSolver<T, Order>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}

#ifdef PHYSICA_MKL
    #include "SymmEigenSolver_MKL.h"
#endif
