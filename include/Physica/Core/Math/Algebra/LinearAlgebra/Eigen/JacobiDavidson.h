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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Orthogonalize.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations/IterateSolver.h"
#include "EigenSolver.h"

namespace Physica::Core {
    /**
     * References:
     * [1] GAMM-Mitteilungen 29(2), 368-382 (2006); https://doi.org/10.1002/gamm.201490038
     * [2] SIAM Review 42(2), 267–293 (2000); https://doi.org/10.1137/S0036144599363084
     * [3] Numerical Linear Algebra with Applications 9(1), 21-44 (2001); https://doi.org/10.1002/nla.246
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
        constexpr static size_t MaxLinearSolverIteration = 64;
        constexpr static double LinearSolverPrecision = 1E-5;
        constexpr static double DefaultStableThreshold = 0.1; //Refer to [3]
        constexpr static double NearConvergeThreshold = 1E-3;
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
        void swap(JacobiDavidson& __restrict other) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return linearSolver.getOrder(); }
        [[nodiscard]] size_t getNumRequired() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const VectorType& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const DenseMatrix<ScalarType>& getEigenvectors() const noexcept { return eigenvectors; }
        /* Setters */
        void setError(RealType error_);
    private:
        template<class MatrixType>
        size_t initSearchSpace(const RValueMatrix<MatrixType>& source,
                                VectorType& initial,
                                size_t eigenIndex);
        void assembleProjects(size_t numSearchDim);
        void ordinarySearch(size_t numSearchDim);
        void refinedSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal);
        template<class MatrixType>
        void prepareRestart(const RValueMatrix<MatrixType>& source, size_t eigenIndex);
        /* Static member */
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
            , error(error_)
            , stableThreshold(stableThreshold_) {
        resize(size, numRequired);
    }

    template<class ScalarType>
    JacobiDavidson<ScalarType>& JacobiDavidson<ScalarType>::operator=(JacobiDavidson<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    template<class MatrixType>
    void JacobiDavidson<ScalarType>::compute(const RValueMatrix<MatrixType>& source, VectorType initial, ScalarType eigenGoal) {
        constexpr bool isHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        constexpr bool isRealSymm = !ScalarType::isComplex && MatrixOption::isSymmMatrix<MatrixType>();
        static_assert(isHermite || isRealSymm, "[Error]: Support for complex eigen problems is not implemented");
        assert(source.getRow() == source.getColumn() && "[Error]: Matrix should be square");
        assert(source.getRow() == initial.getLength() && "[Error]: Dimensions do not match");

        VectorType residule(initial.getLength());
        VectorType buffer(initial.getLength());
        RealType squaredRes;

        for (size_t i = 0; i < getNumRequired(); ++i) {
            const bool useDefaultGoal = i == 0 && eigenGoal == ScalarType(InvalidGoal);
            RealType lastDeltaEigen = std::numeric_limits<ScalarType>::max();
            size_t numSearchDim = initSearchSpace(source, initial, i);

            ScalarType& eigenvalue = eigenvalues[i];
            auto eigenvector = eigenvectors.col(i);
            bool isConverged = false;
            size_t iteration = 0;
        restart:
            while(true) {
                const bool startFromInitial = numSearchDim == 1;
                if (startFromInitial) {
                    eigenvector = searchSpace.col(0).asVector();
                    eigenvalue = projectSearchSpace(0, 0);
                    residule.swap(searchSpace.asArray()[1]);
                    squaredRes = residule.squaredNorm();
                }
                else {
                    if (i == 0) {
                        ordinarySearch(numSearchDim);
                        eigenvalue = eigenSolver.getEigenvalues()[0];
                        residule = source * eigenvector.asVector() - eigenvalue * eigenvector.asVector();
                    }
                    else {
                        refinedSearch(i, numSearchDim, eigenGoal);
                        residule = source * eigenvector.asVector();
                        eigenvalue = eigenvector.asVector().conjugate() * residule;
                        residule -= eigenvalue * eigenvector.asVector();
                    }
                    squaredRes = residule.squaredNorm();
                    isConverged = squaredRes < error && numSearchDim >= MinSearchDim;
                    const bool shouldRestart = numSearchDim == searchSpace.getColumn();
                    if (isConverged || shouldRestart)
                        break;
                }
                /* Update goal */ {
                    if (useDefaultGoal || startFromInitial)
                        eigenGoal = eigenvalue;
                    else {
                        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
                        const VectorType eigenvector2 = subSearchSpace * eigenSolver.getRawEigenvectors().col(numSearchDim - 2);
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
                /* Correction */ {
                    auto new_direction = searchSpace.col(numSearchDim);
                    new_direction = -residule;
                    new_direction.toUnit();
                    linearSolver.solve_functor([this, i, eigenGoal, &buffer, &source](const VectorType& v, VectorType& dot) {
                        auto orthogonalSpace = eigenvectors.leftCols(i + 1);
                        auto head1 = dot.head(i + 1);
                        head1 = orthogonalSpace.hermite() * v;
                        buffer = v - orthogonalSpace * head1;

                        dot = source.getDerived() * buffer - eigenGoal * buffer;
                        auto head2 = buffer.head(i + 1);
                        head2 = orthogonalSpace.hermite() * dot;
                        dot -= orthogonalSpace * head2;
                    }, new_direction);

                    gramSchmidt(searchSpace.leftCols(numSearchDim), new_direction);
                    auto new_dot = dotSpace.col(numSearchDim);
                    new_dot = source * new_direction;
                    assembleProjects(numSearchDim);
                    numSearchDim += 1;
                }
            }

            ++iteration;
            if (!isConverged) {
                if (iteration == MaxIterationPerEigen)
                    throw BadConvergenceException("Exceed max iteration of JacobiDavidson");
                prepareRestart(source, i);
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
            const bool shouldSwap = i != index_min;
            if (shouldSwap) {
                eigenvalues[i].swap(eigenvalues[index_min]);
                eigenvectors.asArray()[i].swap(eigenvectors.asArray()[index_min]);
            }
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::resize(size_t size, size_t numRequired) {
        assert(numRequired < size);
        assert(size > MaxSearchDim);
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
    void JacobiDavidson<ScalarType>::swap(JacobiDavidson<ScalarType>& __restrict other) noexcept {
        assert(this != &other && "[Error]: Self swap is likely a bug");
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
    /**
     * Refer to [1]
     */
    template<class ScalarType>
    template<class MatrixType>
    size_t JacobiDavidson<ScalarType>::initSearchSpace(
            const RValueMatrix<MatrixType>& source,
            VectorType& initial,
            size_t eigenIndex) {
        if (eigenIndex == 0) {
            /* dim = 0 */ {
                initial.toUnit();
                auto col = searchSpace.col(0);
                col = initial;
                auto dot = dotSpace.col(0);
                dot = source.getDerived() * col.asVector();
                assembleProjects(0);

                auto buffer = searchSpace.col(1);
                buffer = dot.asVector() - projectSearchSpace(0, 0) * col.asVector();
                const RealType squaredRes = buffer.squaredNorm();
                const bool isGoodInitial = squaredRes < RealType(NearConvergeThreshold);
                if (isGoodInitial)
                    return 1;
            }
            for (size_t dim = 1; dim < MinSearchDim; ++dim) {
                auto col = searchSpace.col(dim);
                col = source.getDerived() * initial;
                gramSchmidt(searchSpace.leftCols(dim), col);
                initial = col;

                auto dot = dotSpace.col(dim);
                dot = source.getDerived() * col.asVector();
                assembleProjects(dim);
            }
        }
        else {
            const auto lastEigenvector = eigenvectors.col(eigenIndex - 1);
            /* dim = 0 */ {
                auto col = searchSpace.col(0);
                col -= (lastEigenvector.asVector().conjugate() * col.asVector()) * lastEigenvector.asVector();
                col.toUnit();
                auto dot = dotSpace.col(0);
                dot = source.getDerived() * col.asVector();
                assembleProjects(0);
            }
            for (size_t dim = 1; dim < MinSearchDim; ++dim) {
                auto col = searchSpace.col(dim);
                col -= (lastEigenvector.asVector().conjugate() * col.asVector()) * lastEigenvector.asVector();
                gramSchmidt(searchSpace.leftCols(dim), col);

                auto dot = dotSpace.col(dim);
                dot = source.getDerived() * col.asVector();
                assembleProjects(dim);
            }
        }
        return MinSearchDim;
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::assembleProjects(size_t numSearchDim) {
        const size_t i = numSearchDim;
        auto new_direction = searchSpace.col(i);
        auto new_dot = dotSpace.col(i);
        if (i != 0) {
            /* Fill projectSearchSpace */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = projectSearchSpace.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_direction;
                auto row = projectSearchSpace.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                projectSearchSpace(i, i) = new_dot.asVector().conjugate() * new_direction;
            }
            /* Fill projectDotSpace */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = projectDotSpace.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_dot;
                auto row = projectDotSpace.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                projectDotSpace(i, i) = new_dot.squaredNorm();
            }
        }
        else {
            projectSearchSpace(0, 0) = new_dot.asVector().conjugate() * new_direction;
            projectDotSpace(0, 0) = new_dot.squaredNorm();
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::ordinarySearch(size_t numSearchDim) {
        assert(numSearchDim > 1 && "[Error]: No need to search if dim = 1");
        auto corner = projectSearchSpace.topLeftCorner(numSearchDim);
        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner, true);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(0);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(0);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::refinedSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal) {
        auto corner1 = projectSearchSpace.topLeftCorner(numSearchDim);
        auto corner2 = projectDotSpace.topLeftCorner(numSearchDim);
        auto corner = projectSpace.topLeftCorner(numSearchDim);
        corner = eigenGoal * corner1 + eigenGoal.conjugate() * corner1.hermite();
        corner -= corner2;

        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner, true);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(eigenIndex);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(numSearchDim - 1);
    }

    template<class ScalarType>
    template<class MatrixType>
    void JacobiDavidson<ScalarType>::prepareRestart(const RValueMatrix<MatrixType>& source, size_t eigenIndex) {
        size_t dim = 0;
        for (; dim < MinSearchDim - 1; ++dim) {
            const size_t index = searchSpace.getColumn() - (MinSearchDim - 1) + dim;
            searchSpace.asArray()[dim].swap(searchSpace.asArray()[index]);
            dotSpace.asArray()[dim].swap(dotSpace.asArray()[index]);
            assembleProjects(dim);
        }
        searchSpace.asArray()[dim].swap(eigenvectors.asArray()[eigenIndex]);
        gramSchmidt(searchSpace.leftCols(dim), searchSpace.asArray()[dim]);
        auto dot = dotSpace.col(dim);
        dot = source.getDerived() * searchSpace.col(dim);
        assembleProjects(dim);
    }

    template<class ScalarType>
    size_t JacobiDavidson<ScalarType>::calcInnerIteration(size_t order) {
        size_t ite = MaxSizePerMatrix / order;
        if (ite > MaxSearchDim)
            ite = MaxSearchDim;
        else if (ite < MinSearchDim)
            ite = MinSearchDim;
        return ite;
    }

    template<class ScalarType>
    inline void swap(JacobiDavidson<ScalarType>& __restrict obj1, JacobiDavidson<ScalarType>& __restrict obj2) noexcept {
        obj1.swap(obj2);
    }
}
