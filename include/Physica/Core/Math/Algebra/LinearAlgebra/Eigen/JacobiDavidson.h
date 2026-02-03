/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Orthogonalize.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearSystem/IterateSolver.h"
#include "EigenSolver.h"

namespace Physica {
    /**
     * \class JacobiDavidson provides an iterative solution for computing the lowest eigenvalues and their corresponding eigenvectors.
     * If you are interested in obtaining the highest eigenvalues, simply invert the sign of the matrix.
     * 
     * References:
     * [1] GAMM-Mitteilungen 29(2), 368-382 (2006); https://doi.org/10.1002/gamm.201490038
     * [2] SIAM Review 42(2), 267–293 (2000); https://doi.org/10.1137/S0036144599363084
     * [3] Numerical Linear Algebra with Applications 9(1), 21-44 (2001); https://doi.org/10.1002/nla.246
     * [4] Numerical Methods for Solving Large Scale Eigenvalue Problems; https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf
     */
    template<Scalar T>
    class JacobiDavidson {
        using This = JacobiDavidson<T>;
        using Tr = T::RealType;

        constexpr static bool isComplex = T::isComplex;
        constexpr static size_t MaxIterationPerEigen = 4;
        constexpr static size_t MaxSearchDim = 32; //Refer to [1]
        constexpr static size_t MinSearchDim = 5; //Refer to [1]
        constexpr static size_t MaxBufferSize = 1024UL * 1024UL * 1024UL / sizeof(T);
        constexpr static size_t MaxLinearSolverIteration = 64;
        constexpr static double LinearSolverPrecision = 1E-4;
        constexpr static double DefaultStableThreshold = 1E-5;
        constexpr static double NearConvergeThreshold = 1E-3;
    public:
        constexpr static double InvalidGoal = std::numeric_limits<T>::max();
    private:
        IterateSolver<T> linearSolver;
        EigenSolver<T> eigenSolver;

        size_t numSearchDim;
        DenseMatrix<T> searchSpace;
        DenseMatrix<T> dotSpace;
        DenseMatrix<T> searchSpaceProj;
        DenseMatrix<T> dotSpaceProj;
        DenseMatrix<T> spaceProj;

        VectorND<T> eigenvalues;
        DenseMatrix<T> eigenvectors;

        Tr error = std::numeric_limits<Tr>::epsilon();
        Tr stableThreshold = DefaultStableThreshold;
    public:
        JacobiDavidson();
        JacobiDavidson(size_t size, size_t numRequired);
        JacobiDavidson(const This&) = default;
        JacobiDavidson(This&&) noexcept = default;
        ~JacobiDavidson() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source,
                     VectorND<T> initial,
                     T eigenGoal = T(InvalidGoal));
        void sort();
        void resize(size_t size, size_t numRequired);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return eigenvectors.getRow(); }
        [[nodiscard]] size_t getNumRequired() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const auto& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const auto& getEigenvectors() const noexcept { return eigenvectors; }
        /* Setters */
        void setError(Tr error_) noexcept;
        void setStableThreshold(Tr value) noexcept;
    private:
        void initSearchSpace(const Matrix auto& source, VectorND<T>& initial);
        void projSearchSpace(const Matrix auto& source, size_t eigenIndex);
        void assembleProjects(const Matrix auto& source, size_t dim, bool updateDot);

        void extremeSearch();
        void refinedSearch(size_t eigenIndex, T eigenGoal);

        void correction(size_t eigenIndex, T eigenGoal, const Matrix auto& source, VectorND<T>& residule, VectorND<T>& buffer);
        /* Static member */
        static size_t calcSearchSpaceDim(size_t order);
    };

    template<Scalar T>
    JacobiDavidson<T>::JacobiDavidson() : numSearchDim() {
        linearSolver.mustConverge = false;
        linearSolver.setError(LinearSolverPrecision);
        linearSolver.setMaxIteration(MaxLinearSolverIteration);
    }

    template<Scalar T>
    JacobiDavidson<T>::JacobiDavidson(size_t size, size_t numRequired) : JacobiDavidson() {
        linearSolver.resize(size);
        resize(size, numRequired);
    }

    template<Scalar T>
    void JacobiDavidson<T>::compute(const Matrix auto& source, VectorND<T> initial, T eigenGoal) {
        constexpr bool isHermite = source.isStaticHermite();
        constexpr bool isRealSymm = !T::isComplex && source.isStaticSymm();
        static_assert(isHermite || isRealSymm, "[Error]: Support for complex eigenvalues is not implemented");
        assert(source.getRow() == source.getCol() && "[Error]: Matrix should be square");
        assert(source.getRow() == initial.getLength() && "[Error]: Dimensions do not match");
        assert(eigenGoal == InvalidGoal && "[Error]: Not implemented");

        VectorND<T> residule(initial.getLength());
        VectorND<T> buffer(initial.getLength());
        Tr squaredRes;

        for (size_t i = 0; i < getNumRequired(); ++i) {
            const bool isFirstEigen = i == 0;
            if (isFirstEigen)
                initSearchSpace(source, initial);
            else
                projSearchSpace(source, i);

            Tr lastDeltaEigen = std::numeric_limits<T>::max();
            T& eigenvalue = eigenvalues[i];
            auto eigenvector = eigenvectors.col(i);
            bool isConverged = false;
            size_t iteration = 0;
        restart:
            while(true) {
                const bool isGoodInitial = numSearchDim == 1;
                if (isGoodInitial) {
                    eigenvector = searchSpace.col(0);
                    eigenvalue = searchSpaceProj[0, 0];
                    residule = searchSpace.asArray()[1];
                    squaredRes = residule.squaredNorm();
                    eigenGoal = eigenvalue;
                }
                else {
                    if (isFirstEigen) {
                        extremeSearch();
                        eigenvalue = eigenSolver.getEigenvalues()[0].real();
                        residule = source * eigenvector - eigenvalue * eigenvector;
                        eigenGoal = std::min(eigenGoal.real(), eigenvalue.real());
                    }
                    else {
                        refinedSearch(i, eigenGoal);
                        residule = source * eigenvector;
                        eigenvalue = eigenvector.conjugate() * residule;
                        residule -= eigenvalue * eigenvector;
                    }
                    squaredRes = residule.squaredNorm();
                    isConverged = squaredRes < error;
                    const bool shouldRestart = numSearchDim == searchSpace.getCol();
                    if (isConverged || shouldRestart)
                        break;

                    if (!isFirstEigen) { // Update goal
                        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
                        const VectorND<T> eigenvector2 = subSearchSpace * eigenSolver.getRawEigenvectors().col(1);
                        buffer = source * eigenvector2;
                        const T eigenvalue2 = eigenvector2.conjugate() * buffer;
                        const Tr deltaEigen = abs(eigenvalue2.real() - eigenvalue.real());
                        const bool nearConverge = squaredRes < square(deltaEigen);
                        const bool deltaEigenStable = abs(deltaEigen - lastDeltaEigen) < stableThreshold * lastDeltaEigen;
                        const bool goodDeltaEigen = lastDeltaEigen > std::numeric_limits<T>::epsilon();
                        const bool increaseGoal = !goodDeltaEigen || (deltaEigenStable && nearConverge);
                        const bool decreaseGoal = eigenvalue.real() < eigenGoal.real();
                        if (increaseGoal || decreaseGoal)
                            eigenGoal = eigenvalue;
                        lastDeltaEigen = deltaEigen;
                    }
                }
                correction(i, eigenGoal, source, residule, buffer);
            }

            if (!isConverged) {
                iteration += 1;
                if (iteration == MaxIterationPerEigen)
                    throw BadConvergenceException("Exceed max iteration of JacobiDavidson");
                initial = eigenvector;
                residule = source * eigenvector;
                eigenvalue = eigenvector.conjugate() * residule;
                eigenGoal = eigenvalue;
                initSearchSpace(source, initial);
                goto restart;
            }
        }
    }

    template<Scalar T>
    void JacobiDavidson<T>::sort() {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (eigenvalues[j].real() < eigenvalues[index_min].real())
                    index_min = j;
            }

            const bool shouldSwap = i != index_min;
            if (!shouldSwap)
                continue;

            eigenvalues[i].swap(eigenvalues[index_min]);
            eigenvectors.swap_col(i, index_min);
        }
    }

    template<Scalar T>
    void JacobiDavidson<T>::resize(size_t size, size_t numRequired) {
        assert(numRequired <= size && "[Error]: Requiring more eigen pairs than matrix size");
        assert(size > MaxSearchDim && "[Error]: Use dense eigen solver is prefered");
        const size_t dim = calcSearchSpaceDim(size);
        linearSolver.resize(size);
        eigenvalues.resize(numRequired);
        eigenvectors.resize(size, numRequired);
        eigenSolver.resize(dim, true);
        searchSpace.resize(size, dim);
        dotSpace.resize(size, dim);
        searchSpaceProj.resize(dim, dim);
        dotSpaceProj.resize(dim, dim);
        spaceProj.resize(dim, dim);
    }

    template<Scalar T>
    void JacobiDavidson<T>::swap(JacobiDavidson<T>& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        linearSolver.swap(obj.linearSolver);
        eigenSolver.swap(obj.eigenSolver);

        std::swap(numSearchDim, obj.numSearchDim);
        searchSpace.swap(obj.searchSpace);
        dotSpace.swap(obj.dotSpace);
        searchSpaceProj.swap(obj.searchSpaceProj);
        dotSpaceProj.swap(obj.dotSpaceProj);
        spaceProj.swap(obj.spaceProj);

        eigenvalues.swap(obj.eigenvalues);
        eigenvectors.swap(obj.eigenvectors);

        error.swap(obj.error);
        stableThreshold.swap(obj.stableThreshold);
    }

    template<Scalar T>
    void JacobiDavidson<T>::setError(Tr error_) noexcept {
        assert(error_.isPositive() && "[Error]: Invalid argument");
        error = std::move(error_);
    }

    template<Scalar T>
    void JacobiDavidson<T>::setStableThreshold(Tr value) noexcept {
        assert(value.isPositive() && "[Error]: Invalid argument");
        stableThreshold = std::move(value);
    }
    /**
     * Refer to [1]
     */
    template<Scalar T>
    void JacobiDavidson<T>::initSearchSpace(const Matrix auto& source, VectorND<T>& initial) {
        /* dim = 0 */ {
            initial.toUnit();
            auto col = searchSpace.col(0);
            col = initial;
            assembleProjects(source, 0, true);

            const auto dot = dotSpace.col(0);
            auto buffer = searchSpace.col(1);
            buffer = dot - searchSpaceProj[0, 0] * col;
            const Tr squaredRes = buffer.squaredNorm();
            const bool isGoodInitial = squaredRes < Tr(NearConvergeThreshold);
            if (isGoodInitial)
                numSearchDim = 1;
        }
        for (size_t dim = 1; dim < MinSearchDim; ++dim) {
            auto col = searchSpace.col(dim);
            col = source * initial;
            normGramSchmidt(searchSpace.leftCols(dim), col, col.squaredNorm());
            initial = col;
            assembleProjects(source, dim, true);
        }
        numSearchDim = MinSearchDim;
    }
    /**
     * Refer to [1]
     */
    template<Scalar T>
    void JacobiDavidson<T>::projSearchSpace(const Matrix auto& source, size_t eigenIndex) {
        const auto lastEigenvector = eigenvectors.col(eigenIndex - 1);
        for (size_t dim = 0; dim < MinSearchDim; ++dim) {
            auto col = searchSpace.col(dim);
            const T dot1 = lastEigenvector.conjugate() * col;
            const Tr normalizer = reciprocal(sqrt(Tr(1) - dot1.squaredNorm()));
            col -= dot1 * lastEigenvector;
            col *= normalizer;
            if (dim == 0) {
                auto dot2 = dotSpace.col(0);
                dot2 = source * col;
            }
            else
                gramSchmidt(searchSpace.leftCols(dim), col);
            assembleProjects(source, dim, true);
        }
        numSearchDim = MinSearchDim;
    }

    template<Scalar T>
    void JacobiDavidson<T>::assembleProjects(const Matrix auto& source, size_t numSearchDim, bool updateDot) {
        const size_t i = numSearchDim;
        const auto new_direction = searchSpace.col(i);
        auto new_dot = dotSpace.col(i);
        if (updateDot)
            new_dot = source * new_direction;
        if (i != 0) {
            /* Fill searchSpaceProj */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = searchSpaceProj.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_direction;
                auto row = searchSpaceProj.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                searchSpaceProj[i, i] = new_dot.conjugate() * new_direction;
            }
            /* Fill dotSpaceProj */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = dotSpaceProj.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_dot;
                auto row = dotSpaceProj.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                dotSpaceProj[i, i] = new_dot.squaredNorm();
            }
        }
        else {
            searchSpaceProj[0, 0] = new_dot.conjugate() * new_direction;
            dotSpaceProj[0, 0] = new_dot.squaredNorm();
        }
    }

    template<Scalar T>
    void JacobiDavidson<T>::extremeSearch() {
        assert(numSearchDim > 1 && "[Error]: No need to search if dim = 1");
        auto corner = searchSpaceProj.topLeftCorner(numSearchDim);
        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(0);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(0);
    }

    template<Scalar T>
    void JacobiDavidson<T>::refinedSearch(size_t eigenIndex, T eigenGoal) {
        assert(numSearchDim > 1 && "[Error]: No need to search if dim = 1");
        assert(eigenGoal.real() != Tr(InvalidGoal) && "[Error]: Goal is not updated, this is a bug");
        const auto corner1 = dotSpaceProj.topLeftCorner(numSearchDim);
        const auto corner2 = searchSpaceProj.topLeftCorner(numSearchDim);
        auto corner = spaceProj.topLeftCorner(numSearchDim);
        corner = corner1;
        if constexpr (isComplex)
            corner -= eigenGoal * corner2 + eigenGoal.conjugate() * corner2.hermite();
        else
            corner -= T(2) * eigenGoal * corner2;

        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(eigenIndex);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(0);
    }

    template<Scalar T>
    void JacobiDavidson<T>::correction(size_t eigenIndex, T eigenGoal, const Matrix auto& source, VectorND<T>& residule, VectorND<T>& buffer) {
        const auto orthogonalSpace = eigenvectors.leftCols(eigenIndex + 1);
        auto new_direction = searchSpace.col(numSearchDim);
        new_direction = residule;
        gramSchmidt(orthogonalSpace, new_direction);
        linearSolver.cg_functor([this, eigenIndex, eigenGoal, &buffer, &source](const VectorND<T>& v, VectorND<T>& dot) {
            const auto orthogonalSpace = eigenvectors.leftCols(eigenIndex + 1);
            auto head1 = dot.head(eigenIndex + 1);
            head1 = orthogonalSpace.hermite() * v;
            buffer = v - orthogonalSpace * head1;

            dot = source * buffer - eigenGoal * buffer;
            auto head2 = buffer.head(eigenIndex + 1);
            head2 = orthogonalSpace.hermite() * dot;
            dot -= orthogonalSpace * head2;
        }, new_direction);

        new_direction.toUnit();
        gramSchmidt(searchSpace.leftCols(numSearchDim), new_direction);
        assembleProjects(source, numSearchDim, true);
        numSearchDim += 1;
    }

    template<Scalar T>
    size_t JacobiDavidson<T>::calcSearchSpaceDim(size_t order) {
        const size_t result = MaxBufferSize / order;
        if (result > MaxSearchDim)
            return MaxSearchDim;
        if (result < MinSearchDim)
            return MinSearchDim;
        return result;
    }
}

namespace std {
    template<Physica::Scalar T>
    void swap(
            Physica::JacobiDavidson<T>& __restrict obj1,
            Physica::JacobiDavidson<T>& __restrict obj2) noexcept {
        obj1.swap(obj2);
    }
}
