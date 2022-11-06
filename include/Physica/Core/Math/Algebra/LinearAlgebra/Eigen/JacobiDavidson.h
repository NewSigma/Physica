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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations/IterateSolver.h"
#include "EigenSolver.h"

namespace Physica::Core {
    /**
     * References:
     * [1] M.E. Hochstenbach and Y. Notay, The Jacobi–Davidson method (https://doi.org/10.1002/gamm.201490038)
     * [2] Gerard L. G. Sleijpen and Henk A. Van der Vorst, A Jacobi-Davidson Iteration Method for Linear Eigenvalue Problems (https://doi.org/10.1137/S0036144599363084)
     * [3] https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf
     * [4] Y. Notay, Combination of Jacobi–Davidson and conjugate gradients for the partial symmetric eigenproblem (https://doi.org/10.1002/nla.246)
     */
    template<class ScalarType>
    class JacobiDavidson {
        using RealType = typename ScalarType::RealType;
        using WorkingMatrix = DenseMatrix<ScalarType>;
        using VectorType = Vector<ScalarType>;
        using LinearSolverType = IterateSolver<ScalarType>;
        using EigenSolverType = EigenSolver<DenseMatrix<ScalarType>>;

        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static size_t MaxIterationPerEigen = 4;
        constexpr static size_t MaxSearchDim = 32; //Refer to [1]
        constexpr static size_t MinSearchDim = 5; //Refer to [1]
        constexpr static size_t MaxSizePerMatrix = 1024 * 1024 * 1024 / sizeof(ScalarType);
        constexpr static size_t MaxLinearSolverIteration = 0;
        constexpr static double LinearSolverPrecision = 1E-5;
        constexpr static double DefaultStableThreshold = 0.1; //Refer to [4]
    public:
        constexpr static double InvalidGoal = std::numeric_limits<ScalarType>::max();
    private:
        LinearSolverType linearSolver;
        EigenSolverType eigenSolver;
        DenseMatrix<ScalarType> searchSpace;
        DenseMatrix<ScalarType> dotSpace;
        DenseMatrix<ScalarType> projectSearchSpace;
        DenseMatrix<ScalarType> projectDotSpace;
        DenseMatrix<ScalarType> projectSpace;
        VectorType eigenvalues;
        DenseMatrix<ScalarType> eigenvectors;
        RealType error;
        RealType stableThreshold;
    public:
        JacobiDavidson();
        JacobiDavidson(size_t size,
                       size_t numRequired,
                       RealType error_ = RealType(std::numeric_limits<RealType>::epsilon()),
                       RealType stableThrehold_ = DefaultStableThreshold);
        JacobiDavidson(const JacobiDavidson&) = default;
        JacobiDavidson(JacobiDavidson&&) noexcept = default;
        ~JacobiDavidson() = default;
        /* Operators */
        JacobiDavidson& operator=(JacobiDavidson obj) noexcept;
        /* Operations */
        template<class MatrixType>
        void compute(const RValueMatrix<MatrixType>& source,
                     VectorType initial,
                     ScalarType eigenGoal = ScalarType(InvalidGoal));
        void sort();
        void resize(size_t size, size_t numRequired);
        void swap(JacobiDavidson& other) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return linearSolver.getOrder(); }
        [[nodiscard]] size_t getNumRequired() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const VectorType& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const DenseMatrix<ScalarType>& getEigenvectors() const noexcept { return eigenvectors; }
        /* Setters */
        void setError(RealType error_);
    private:
        void assembleProjects(size_t numSearchDim);
        void refinedRitzSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal);
        static size_t calcInnerIteration(size_t order);
    };

    template<class ScalarType>
    JacobiDavidson<ScalarType>::JacobiDavidson()
            : linearSolver(MaxLinearSolverIteration, LinearSolverPrecision)
            , error(std::numeric_limits<RealType>::epsilon())
            , stableThreshold(DefaultStableThreshold) {}

    template<class ScalarType>
    JacobiDavidson<ScalarType>::JacobiDavidson(size_t size, size_t numRequired, RealType error_, RealType stableThreshold_)
            : linearSolver(MaxLinearSolverIteration, LinearSolverPrecision)
            , eigenvalues(numRequired)
            , eigenvectors(size, numRequired)
            , error(error_)
            , stableThreshold(stableThreshold_) {
        assert(numRequired < size);
        const size_t ite = calcInnerIteration(size);
        eigenSolver.resize(ite);
        searchSpace.resize(size, ite);
        dotSpace.resize(size, ite);
        projectSearchSpace.resize(ite, ite);
        projectDotSpace.resize(ite, ite);
        projectSpace.resize(ite, ite);
    }

    template<class ScalarType>
    JacobiDavidson<ScalarType>& JacobiDavidson<ScalarType>::operator=(JacobiDavidson<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    template<class MatrixType>
    void JacobiDavidson<ScalarType>::compute(const RValueMatrix<MatrixType>& source, VectorType initial, ScalarType eigenGoal) {
        static_assert(std::is_same<MatrixType, DenseHermiteMatrix<ScalarType>>::value || std::is_same<MatrixType, DenseSymmMatrix<ScalarType>>::value, "[Error]: Not implemented");
        assert(source.getRow() == source.getColumn());
        assert(source.getRow() == initial.getLength());

        VectorType residule;
        VectorType buffer(initial.getLength());

        for (size_t i = 0; i < getNumRequired(); ++i) {
            const bool useDefaultGoal = i == 0 && eigenGoal == ScalarType(InvalidGoal);
            RealType lastDeltaEigen = std::numeric_limits<ScalarType>::max();
            ScalarType& eigenvalue = eigenvalues[i];
            auto eigenvector = eigenvectors.col(i);
            /* Init */ {
                if (i == 0)
                    eigenvector = initial;
                else {
                    buffer = source.getDerived() * initial;
                    buffer.swap(initial);
                    initial.toUnit();
                    auto orthogonalSpace = eigenvectors.leftCols(i);
                    auto head = buffer.head(i);
                    head = orthogonalSpace.transpose().conjugate() * initial;
                    eigenvector = initial - orthogonalSpace * head;
                }
                eigenvector.toUnit();

                auto col = searchSpace.col(0);
                col = eigenvector.asVector();
                auto dot = dotSpace.col(0);
                dot = source * eigenvector.asVector();
                eigenvalue = eigenvector.asVector().conjugate() * dot;
                residule = dot - eigenvalue * eigenvector.asVector();
                projectSearchSpace(0, 0) = eigenvalue;
                projectDotSpace(0, 0) = dot.squaredNorm();

                if (useDefaultGoal)
                    eigenGoal = eigenvalue;
            }

            bool isConverged = false;
            size_t numSearchDim = 1;
            size_t iteration = 0;
        restart:
            for (; numSearchDim < searchSpace.getColumn(); ++numSearchDim) {
                /* Correction */ {
                    residule = -residule;
                    residule.toUnit();
                    linearSolver.solve_functor([this, i, eigenGoal, &buffer, &source](const VectorType& v, VectorType& dot) {
                        auto orthogonalSpace = eigenvectors.leftCols(i + 1);
                        auto head1 = dot.head(i + 1);
                        head1 = orthogonalSpace.transpose().conjugate() * v;
                        buffer = v - orthogonalSpace * head1;

                        dot = source.getDerived() * buffer - eigenGoal * buffer;
                        auto head2 = buffer.head(i + 1);
                        head2 = orthogonalSpace.transpose().conjugate() * dot;
                        dot -= orthogonalSpace * head2;
                    }, residule);
                }
                /* Orthogonalize */ {
                    const VectorType& correction = residule;
                    auto leftCols = searchSpace.leftCols(numSearchDim);
                    const VectorType project = leftCols.transpose().conjugate() * correction;
                    auto new_direction = searchSpace.col(numSearchDim);
                    new_direction = correction - leftCols * project;
                    new_direction.toUnit();
                    auto new_dot = dotSpace.col(numSearchDim);
                    new_dot = source * new_direction;
                }
                refinedRitzSearch(i, numSearchDim, eigenGoal);
                residule = source * eigenvector.asVector();
                eigenvalue = eigenvector.asVector().conjugate() * residule;
                residule -= eigenvalue * eigenvector.asVector();

                const RealType squaredRes = residule.squaredNorm();
                isConverged = squaredRes < error;
                if (isConverged)
                    break;
                /* Update goal */ {
                    if (useDefaultGoal)
                        eigenGoal = eigenvalue;
                    else {
                        auto subSearchSpace = searchSpace.leftCols(numSearchDim + 1);
                        const VectorType eigenvector2 = subSearchSpace * eigenSolver.getRawEigenvectors().col(numSearchDim - 1);
                        buffer = source * eigenvector2;
                        const ScalarType eigenvalue2 = eigenvector2.conjugate() * buffer;
                        const RealType deltaEigen = eigenvalue2.getReal() - eigenvalue.getReal();
                        const bool nearConverge = squaredRes < square(deltaEigen);
                        const bool deltaEigenStable = abs(deltaEigen / lastDeltaEigen - 1.0) < stableThreshold;
                        const bool goodDeltaEigen = lastDeltaEigen > std::numeric_limits<ScalarType>::epsilon();
                        if (!goodDeltaEigen || (deltaEigenStable && nearConverge))
                            eigenGoal = eigenvalue;
                        lastDeltaEigen = deltaEigen;
                    }
                }
            }

            ++iteration;
            if (!isConverged) {
                if (iteration == MaxIterationPerEigen)
                    throw BadConvergenceException("Exceed max iteration of JacobiDavidson");
                numSearchDim = MinSearchDim;
                goto restart;
            }
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::sort() {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (eigenvalues[j].getReal() < eigenvalues[index_min].getReal())
                    index_min = j;
            }
            eigenvalues[i].swap(eigenvalues[index_min]);
            eigenvectors[i].swap(eigenvectors[index_min]);
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::resize(size_t size, size_t numRequired) {
        const size_t ite = calcInnerIteration(size);
        eigenvalues.resize(numRequired);
        eigenvectors.resize(size, numRequired);
        eigenSolver.resize(ite);
        searchSpace.resize(size, ite);
        dotSpace.resize(size, ite);
        projectSearchSpace.resize(ite, ite);
        projectDotSpace.resize(ite, ite);
        projectSpace.resize(ite, ite);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::swap(JacobiDavidson<ScalarType>& other) noexcept {
        linearSolver.swap(other.linearSolver);
        eigenSolver.swap(other.eigenSolver);
        searchSpace.swap(other.searchSpace);
        dotSpace.swap(other.dotSpace);
        projectSearchSpace.swap(other.projectSearchSpace);
        projectDotSpace.swap(other.projectDotSpace);
        projectSpace.swap(other.projectSpace);
        eigenvalues.swap(other.eigenvalues);
        eigenvectors.swap(other.eigenvectors);
        error.swap(other.error);
        stableThreshold.swap(other.stableThreshold);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::setError(RealType error_) {
        assert(error_.isPositive());
        error = error_;
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::assembleProjects(size_t numSearchDim) {
        const size_t i = numSearchDim;
        auto new_direction = searchSpace.col(i);
        auto new_dot = dotSpace.col(i);
        /* Fill projectSearchSpace */ {
            auto leftCols = dotSpace.leftCols(i);
            auto col = projectSearchSpace.col(i);
            auto head1 = col.head(i);
            head1 = leftCols.transpose().conjugate() * new_direction;
            auto row = projectSearchSpace.row(i);
            auto head2 = row.head(i);
            head2 = head1.conjugate();
            projectSearchSpace(i, i) = new_dot.asVector().conjugate() * new_direction;
        }
        /* Fill projectDotSpace */ {
            auto leftCols = dotSpace.leftCols(i);
            auto col = projectDotSpace.col(i);
            auto head1 = col.head(i);
            head1 = leftCols.transpose().conjugate() * new_dot;
            auto row = projectDotSpace.row(i);
            auto head2 = row.head(i);
            head2 = head1.conjugate();
            projectDotSpace(i, i) = new_dot.squaredNorm();
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::refinedRitzSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal) {
        const size_t i = numSearchDim;
        assembleProjects(i);

        auto corner1 = projectSearchSpace.topLeftCorner(i + 1);
        auto corner2 = projectDotSpace.topLeftCorner(i + 1);
        auto corner = projectSpace.topLeftCorner(i + 1);
        corner = eigenGoal * corner1 + eigenGoal.conjugate() * corner1.transpose().conjugate();
        corner -= corner2;

        eigenSolver.resize(i + 1);
        eigenSolver.compute(corner, true);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(i + 1);
        auto eigenvector = eigenvectors.col(eigenIndex);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(i);
    }

    template<class ScalarType>
    size_t JacobiDavidson<ScalarType>::calcInnerIteration(size_t order) {
        size_t ite = MaxSizePerMatrix / order;
        if (ite > MaxSearchDim)
            ite = MaxSearchDim;
        else if (ite < MinSearchDim)
            ite = MinSearchDim;
        return std::min(ite, order);
    }

    template<class ScalarType>
    inline void swap(JacobiDavidson<ScalarType>& obj1, JacobiDavidson<ScalarType>& obj2) noexcept {
        obj1.swap(obj2);
    }
}
