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
     */
    template<class ScalarType>
    class JacobiDavidson {
        using RealType = typename ScalarType::RealType;
        using WorkingMatrix = DenseMatrix<ScalarType>;
        using VectorType = Vector<ScalarType>;
        using LinearSolverType = IterateSolver<ScalarType>;
        using EigenSolverType = EigenSolver<DenseMatrix<ScalarType>>;

        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static size_t MaxInnerIteration = 32; //Refer to [1]
        constexpr static size_t MinInnerIteration = 4; //Refer to [1]
        constexpr static size_t MaxSizePerMatrix = 1024 * 1024 * 1024 / sizeof(ScalarType);
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
    public:
        JacobiDavidson() = default;
        JacobiDavidson(size_t size, size_t numRequired, RealType error_ = RealType(std::numeric_limits<RealType>::epsilon()));
        JacobiDavidson(const JacobiDavidson&) = default;
        JacobiDavidson(JacobiDavidson&&) noexcept = default;
        ~JacobiDavidson() = default;
        /* Operators */
        JacobiDavidson& operator=(JacobiDavidson obj) noexcept;
        /* Operations */
        template<class MatrixType>
        void compute(const RValueMatrix<MatrixType>& source, VectorType initial);
        void swap(JacobiDavidson& other) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return linearSolver.getOrder(); }
        [[nodiscard]] size_t getInnerIteration() const noexcept { return searchSpace.getColumn(); }
        [[nodiscard]] size_t getNumRequired() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const VectorType& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const DenseMatrix<ScalarType>& getEigenvectors() const noexcept { return eigenvectors; }
    private:
        void refinedRitzSearch(size_t index, size_t iteration);
        static size_t calcInnerIteration(size_t order);
    };

    template<class ScalarType>
    JacobiDavidson<ScalarType>::JacobiDavidson(size_t size, size_t numRequired, RealType error_)
            : eigenvalues(numRequired), eigenvectors(size, numRequired), error(error_) {
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
    void JacobiDavidson<ScalarType>::compute(const RValueMatrix<MatrixType>& source, VectorType initial) {
        static_assert(std::is_same<MatrixType, DenseHermiteMatrix<ScalarType>>::value || std::is_same<MatrixType, DenseSymmMatrix<ScalarType>>::value, "[Error]: Not implemented");
        assert(source.getRow() == source.getColumn());
        assert(source.getRow() == initial.getLength());

        VectorType residule;
        VectorType buffer(initial.getLength());

        const size_t innerIteration = getInnerIteration();
        for (size_t i = 0; i < getNumRequired(); ++i) {
            ScalarType& eigenvalue = eigenvalues[i];
            auto eigenvector = eigenvectors.col(i);
            if (i == 0)
                eigenvector = initial;
            else {
                buffer = source.getDerived() * initial;
                buffer.swap(initial);
                eigenvector = initial;
                auto orthogonalSpace = eigenvectors.leftCols(i);
                auto head = buffer.head(i);
                head = orthogonalSpace.transpose().conjugate() * eigenvector;
                eigenvector -= orthogonalSpace * head;
            }
            eigenvector.toUnit();
            bool isConverged = false;
        restart:
            /* j = 0 */ {
                auto col = searchSpace.col(0);
                col = eigenvector.asVector();
                auto dot = dotSpace.col(0);
                dot = source * eigenvector.asVector();
                eigenvalue = eigenvector.asVector().conjugate() * dot;
                residule = dot - eigenvalue * eigenvector.asVector();
                projectSearchSpace(0, 0) = eigenvalue;
                projectDotSpace(0, 0) = dot.squaredNorm();
            }

            for (size_t j = 1; j < innerIteration; ++j) {
                /* Correction */ {
                    residule = -residule;
                    linearSolver.solve_functor([this, i, &buffer, &source](const VectorType& v, VectorType& dot) {
                        auto orthogonalSpace = eigenvectors.leftCols(i + 1);
                        auto head1 = dot.head(i + 1);
                        head1 = orthogonalSpace.transpose().conjugate() * v;
                        buffer = v - orthogonalSpace * head1;

                        dot = source.getDerived() * buffer - eigenvalues[i] * buffer;
                        auto head2 = buffer.head(i + 1);
                        head2 = orthogonalSpace.transpose().conjugate() * dot;
                        dot -= orthogonalSpace * head2;
                    }, residule);
                }
                /* Orthogonalize */ {
                    const VectorType& correction = residule;
                    auto leftCols = searchSpace.leftCols(j);
                    const VectorType project = leftCols.transpose().conjugate() * correction;
                    auto new_direction = searchSpace.col(j);
                    new_direction = correction - leftCols * project;
                    new_direction.toUnit();
                    auto new_dot = dotSpace.col(j);
                    new_dot = source * new_direction;
                }
                refinedRitzSearch(i, j);
                residule = source * eigenvector.asVector();
                eigenvalue = eigenvector.asVector().conjugate() * residule;
                residule -= eigenvalue * eigenvector.asVector();

                isConverged = residule.squaredNorm() < error;
                if (isConverged)
                    break;
            }

            if (!isConverged)
                goto restart;
        }
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
        error.swap(other.error);
        eigenvalues.swap(other.eigenvalues);
        eigenvectors.swap(other.eigenvectors);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::refinedRitzSearch(size_t index, size_t iteration) {
        const size_t j = iteration;
        auto new_direction = searchSpace.col(j);
        auto new_dot = dotSpace.col(j);
        /* Fill projectSearchSpace */ {
            auto leftCols = dotSpace.leftCols(j);
            auto col = projectSearchSpace.col(j);
            auto head1 = col.head(j);
            head1 = leftCols.transpose().conjugate() * new_direction;
            auto row = projectSearchSpace.row(j);
            auto head2 = row.head(j);
            head2 = head1.conjugate();
            projectSearchSpace(j, j) = new_dot.asVector().conjugate() * new_direction;
        }
        /* Fill projectDotSpace */ {
            auto leftCols = dotSpace.leftCols(j);
            auto col = projectDotSpace.col(j);
            auto head1 = col.head(j);
            head1 = leftCols.transpose().conjugate() * new_dot;
            auto row = projectDotSpace.row(j);
            auto head2 = row.head(j);
            head2 = head1.conjugate();
            projectDotSpace(j, j) = new_dot.squaredNorm();
        }
        /* Search */ {
            auto corner1 = projectSearchSpace.topLeftCorner(j + 1);
            auto corner2 = projectDotSpace.topLeftCorner(j + 1);
            auto corner = projectSpace.topLeftCorner(j + 1);

            const ScalarType target = toRealVector(eigenvalues.head(index + 1)).min();
            corner = target * corner1 + target.conjugate() * corner1.transpose().conjugate();
            corner -= corner2;

            eigenSolver.resize(j + 1);
            eigenSolver.compute(corner, true);
            eigenSolver.sort();
            auto subSearchSpace = searchSpace.leftCols(j + 1);
            auto eigenvector = eigenvectors.col(index);
            eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(j);
        }
    }

    template<class ScalarType>
    size_t JacobiDavidson<ScalarType>::calcInnerIteration(size_t order) {
        size_t ite = MaxSizePerMatrix / order;
        if (ite > MaxInnerIteration)
            ite = MaxInnerIteration;
        else if (ite < MinInnerIteration)
            ite = MinInnerIteration;
        return std::min(ite, order);
    }

    template<class ScalarType>
    inline void swap(JacobiDavidson<ScalarType>& obj1, JacobiDavidson<ScalarType>& obj2) noexcept {
        obj1.swap(obj2);
    }
}
