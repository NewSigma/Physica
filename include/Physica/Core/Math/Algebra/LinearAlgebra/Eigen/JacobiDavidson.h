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
     * \class JacobiDavidson provides an iterative solution for computing the lowest eigenvalues and their corresponding eigenvectors.
     * If you are interested in obtaining the highest eigenvalues, simply invert the sign of the matrix.
     * 
     * References:
     * [1] GAMM-Mitteilungen 29(2), 368-382 (2006); https://doi.org/10.1002/gamm.201490038
     * [2] SIAM Review 42(2), 267–293 (2000); https://doi.org/10.1137/S0036144599363084
     * [3] Numerical Linear Algebra with Applications 9(1), 21-44 (2001); https://doi.org/10.1002/nla.246
     * [4] Numerical Methods for Solving Large Scale Eigenvalue Problems; https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf
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
        constexpr static size_t MaxBufferSize = 1024 * 1024 * 1024 / sizeof(ScalarType);
        constexpr static size_t MaxLinearSolverIteration = 64;
        constexpr static double LinearSolverPrecision = 1E-4;
        constexpr static double DefaultStableThreshold = 1E-5;
        constexpr static double NearConvergeThreshold = 1E-3;
    public:
        constexpr static double InvalidGoal = std::numeric_limits<ScalarType>::max();
    private:
        LinearSolverType linearSolver;
        EigenSolverType eigenSolver;
        DenseMatrix<ScalarType> searchSpace;
        DenseMatrix<ScalarType> dotSpace;
        DenseMatrix<ScalarType> searchSpaceProj;
        DenseMatrix<ScalarType> dotSpaceProj;
        DenseMatrix<ScalarType> spaceProj;
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
        void compute(const RValueMatrix<MatrixType>& source_,
                     VectorType initial,
                     ScalarType eigenGoal = ScalarType(InvalidGoal));
        void sort();
        void resize(size_t size, size_t numRequired);
        void swap(JacobiDavidson& __restrict other) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return eigenvectors.getRow(); }
        [[nodiscard]] size_t getNumRequired() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const VectorType& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const DenseMatrix<ScalarType>& getEigenvectors() const noexcept { return eigenvectors; }
        /* Setters */
        inline void setError(RealType error_) noexcept;
    private:
        template<class MatrixType>
        size_t initSearchSpace(const RValueMatrix<MatrixType>& source_, VectorType& initial);
        template<class MatrixType>
        size_t projSearchSpace(const RValueMatrix<MatrixType>& source_, size_t eigenIndex);
        template<class MatrixType>
        void assembleProjects(const RValueMatrix<MatrixType>& source_, size_t numSearchDim, bool updateDot);
        void extremeSearch(size_t numSearchDim);
        void refinedSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal);
        /* Static member */
        static size_t calcSearchSpaceDim(size_t order);
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
    void JacobiDavidson<ScalarType>::compute(const RValueMatrix<MatrixType>& source_, VectorType initial, ScalarType eigenGoal) {
        constexpr bool isHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        constexpr bool isRealSymm = !ScalarType::isComplex && MatrixOption::isSymmMatrix<MatrixType>();
        static_assert(isHermite || isRealSymm, "[Error]: Support for complex eigenvalues is not implemented");
        const auto& source = source_.getDerived();
        assert(source.getRow() == source.getColumn() && "[Error]: Matrix should be square");
        assert(source.getRow() == initial.getLength() && "[Error]: Dimensions do not match");

        VectorType residule(initial.getLength());
        VectorType buffer(initial.getLength());
        RealType squaredRes;

        for (size_t i = 0; i < getNumRequired(); ++i) {
            const bool useDefaultGoal = i == 0 && eigenGoal == ScalarType(InvalidGoal);
            RealType lastDeltaEigen = std::numeric_limits<ScalarType>::max();
            size_t numSearchDim = i == 0 ? initSearchSpace(source, initial) : projSearchSpace(source, i);

            ScalarType& eigenvalue = eigenvalues[i];
            auto eigenvector = eigenvectors.col(i);
            bool isConverged = false;
            size_t iteration = 0;
        restart:
            while(true) {
                const bool startFromInitial = numSearchDim == 1;
                if (startFromInitial) {
                    eigenvector = searchSpace.col(0).asVector();
                    eigenvalue = searchSpaceProj(0, 0);
                    residule.swap(searchSpace.asArray()[1]);
                    squaredRes = residule.squaredNorm();
                }
                else {
                    if (i == 0) {
                        extremeSearch(numSearchDim);
                        eigenvalue = eigenSolver.getEigenvalues()[0].getReal();
                        residule = source * eigenvector.asVector() - eigenvalue * eigenvector.asVector();
                    }
                    else {
                        refinedSearch(i, numSearchDim, eigenGoal);
                        residule = source * eigenvector.asVector();
                        eigenvalue = eigenvector.asVector().conjugate() * residule;
                        residule -= eigenvalue * eigenvector.asVector();
                    }
                    squaredRes = residule.squaredNorm();
                    isConverged = squaredRes < error;
                    const bool shouldRestart = numSearchDim == searchSpace.getColumn();
                    if (isConverged || shouldRestart)
                        break;
                }
                /* Update goal */ {
                    if (useDefaultGoal || startFromInitial)
                        eigenGoal = eigenvalue;
                    else {
                        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
                        const VectorType eigenvector2 = subSearchSpace * eigenSolver.getRawEigenvectors().col(1);
                        buffer = source * eigenvector2;
                        const ScalarType eigenvalue2 = eigenvector2.conjugate() * buffer;
                        const RealType deltaEigen = abs(eigenvalue2.getReal() - eigenvalue.getReal());
                        const bool nearConverge = squaredRes < square(deltaEigen);
                        const bool deltaEigenStable = abs(deltaEigen - lastDeltaEigen) < stableThreshold * lastDeltaEigen;
                        const bool goodDeltaEigen = lastDeltaEigen > std::numeric_limits<ScalarType>::epsilon();
                        const bool increaseGoal = !goodDeltaEigen || (deltaEigenStable && nearConverge);
                        const bool decreaseGoal = eigenvalue.getReal() < eigenGoal.getReal();
                        if (increaseGoal || decreaseGoal)
                            eigenGoal = eigenvalue;
                        lastDeltaEigen = deltaEigen;
                    }
                }
                /* Correction */ {
                    const auto orthogonalSpace = eigenvectors.leftCols(i + 1);
                    auto new_direction = searchSpace.col(numSearchDim);
                    new_direction = residule;
                    normGramSchmidt(orthogonalSpace, new_direction, squaredRes);
                    linearSolver.solve_functor([this, i, eigenGoal, &buffer, &source](const VectorType& v, VectorType& dot) {
                        const auto orthogonalSpace = eigenvectors.leftCols(i + 1);
                        auto head1 = dot.head(i + 1);
                        head1 = orthogonalSpace.hermite() * v;
                        buffer = v - orthogonalSpace * head1;

                        dot = source * buffer - eigenGoal * buffer;
                        auto head2 = buffer.head(i + 1);
                        head2 = orthogonalSpace.hermite() * dot;
                        dot -= orthogonalSpace * head2;
                    }, new_direction);

                    new_direction.toUnit();
                    normGramSchmidt(searchSpace.leftCols(numSearchDim), new_direction);
                    assembleProjects(source, numSearchDim, true);
                    numSearchDim += 1;
                }
            }

            if (!isConverged) {
                iteration += 1;
                if (iteration == MaxIterationPerEigen)
                    throw BadConvergenceException("Exceed max iteration of JacobiDavidson");
                initial = eigenvector;
                residule = source * eigenvector.asVector();
                eigenvalue = eigenvector.asVector().conjugate() * residule;
                eigenGoal = eigenvalue;
                numSearchDim = initSearchSpace(source_, initial);
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
            if (!shouldSwap)
                continue;

            eigenvalues[i].swap(eigenvalues[index_min]);
            eigenvectors.asArray()[i].swap(eigenvectors.asArray()[index_min]);
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::resize(size_t size, size_t numRequired) {
        assert(numRequired <= size && "[Error]: Requiring more eigen pairs than matrix size");
        assert(size > MaxSearchDim && "[Error]: Use dense eigen solver is prefered");
        const size_t dim = calcSearchSpaceDim(size);
        eigenvalues.resize(numRequired);
        eigenvectors.resize(size, numRequired);
        eigenSolver.resize(dim);
        searchSpace.resize(size, dim);
        dotSpace.resize(size, dim);
        searchSpaceProj.resize(dim, dim);
        dotSpaceProj.resize(dim, dim);
        spaceProj.resize(dim, dim);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::swap(JacobiDavidson<ScalarType>& __restrict other) noexcept {
        assert(this != &other && "[Error]: Self swap is likely a bug");
        linearSolver.swap(other.linearSolver);
        eigenSolver.swap(other.eigenSolver);
        searchSpace.swap(other.searchSpace);
        dotSpace.swap(other.dotSpace);
        searchSpaceProj.swap(other.searchSpaceProj);
        dotSpaceProj.swap(other.dotSpaceProj);
        spaceProj.swap(other.spaceProj);
        eigenvalues.swap(other.eigenvalues);
        eigenvectors.swap(other.eigenvectors);
        error.swap(other.error);
        stableThreshold.swap(other.stableThreshold);
    }

    template<class ScalarType>
    inline void JacobiDavidson<ScalarType>::setError(RealType error_) noexcept {
        assert(error_.isPositive() && "[Error]: Invalid argument");
        error = std::move(error_);
    }
    /**
     * Refer to [1]
     */
    template<class ScalarType>
    template<class MatrixType>
    size_t JacobiDavidson<ScalarType>::initSearchSpace(const RValueMatrix<MatrixType>& source_, VectorType& initial) {
        const auto& source = source_.getDerived();
        /* dim = 0 */ {
            initial.toUnit();
            auto col = searchSpace.col(0);
            col = initial;
            assembleProjects(source_, 0, true);

            const auto dot = dotSpace.col(0);
            auto buffer = searchSpace.col(1);
            buffer = dot.asVector() - searchSpaceProj(0, 0) * col.asVector();
            const RealType squaredRes = buffer.squaredNorm();
            const bool isGoodInitial = squaredRes < RealType(NearConvergeThreshold);
            if (isGoodInitial)
                return 1;
        }
        for (size_t dim = 1; dim < MinSearchDim; ++dim) {
            auto col = searchSpace.col(dim);
            col = source * initial;
            normGramSchmidt(searchSpace.leftCols(dim), col, col.squaredNorm());
            initial = col;
            assembleProjects(source_, dim, true);
        }
        return MinSearchDim;
    }
    /**
     * Refer to [1]
     */
    template<class ScalarType>
    template<class MatrixType>
    size_t JacobiDavidson<ScalarType>::projSearchSpace(const RValueMatrix<MatrixType>& source_, size_t eigenIndex) {
        const auto& source = source_.getDerived();
        const auto lastEigenvector = eigenvectors.col(eigenIndex - 1);
        for (size_t dim = 0; dim < MinSearchDim; ++dim) {
            auto col = searchSpace.col(dim);
            const ScalarType dot1 = lastEigenvector.asVector().conjugate() * col.asVector();
            const RealType normalizer = reciprocal(sqrt(RealType(1) - dot1.squaredNorm()));
            col -= dot1 * lastEigenvector.asVector();
            col *= normalizer;
            if (dim == 0) {
                auto dot2 = dotSpace.col(0);
                dot2 = source * col.asVector();
            }
            else
                gramSchmidt(searchSpace.leftCols(dim), col);
            assembleProjects(source_, dim, true);
        }
        return MinSearchDim;
    }

    template<class ScalarType>
    template<class MatrixType>
    void JacobiDavidson<ScalarType>::assembleProjects(const RValueMatrix<MatrixType>& source_, size_t numSearchDim, bool updateDot) {
        const size_t i = numSearchDim;
        const auto new_direction = searchSpace.col(i);
        auto new_dot = dotSpace.col(i);
        if (updateDot)
            new_dot = source_.getDerived() * new_direction.asVector();
        if (i != 0) {
            /* Fill searchSpaceProj */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = searchSpaceProj.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_direction;
                auto row = searchSpaceProj.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                searchSpaceProj(i, i) = new_dot.asVector().conjugate() * new_direction;
            }
            /* Fill dotSpaceProj */ {
                auto leftCols = dotSpace.leftCols(i);
                auto col = dotSpaceProj.col(i);
                auto head1 = col.head(i);
                head1 = leftCols.hermite() * new_dot;
                auto row = dotSpaceProj.row(i);
                auto head2 = row.head(i);
                head2 = head1.conjugate();
                dotSpaceProj(i, i) = new_dot.squaredNorm();
            }
        }
        else {
            searchSpaceProj(0, 0) = new_dot.asVector().conjugate() * new_direction;
            dotSpaceProj(0, 0) = new_dot.squaredNorm();
        }
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::extremeSearch(size_t numSearchDim) {
        assert(numSearchDim > 1 && "[Error]: No need to search if dim = 1");
        auto corner = searchSpaceProj.topLeftCorner(numSearchDim);
        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner, true);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(0);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(0);
    }

    template<class ScalarType>
    void JacobiDavidson<ScalarType>::refinedSearch(size_t eigenIndex, size_t numSearchDim, ScalarType eigenGoal) {
        assert(numSearchDim > 1 && "[Error]: No need to search if dim = 1");
        const auto corner1 = dotSpaceProj.topLeftCorner(numSearchDim);
        const auto corner2 = searchSpaceProj.topLeftCorner(numSearchDim);
        auto corner = spaceProj.topLeftCorner(numSearchDim);
        corner = corner1;
        if constexpr (isComplex)
            corner -= eigenGoal * corner2 + eigenGoal.conjugate() * corner2.hermite();
        else
            corner -= ScalarType(2) * eigenGoal * corner2;

        eigenSolver.resize(numSearchDim);
        eigenSolver.compute(corner, true);
        eigenSolver.sort();
        auto subSearchSpace = searchSpace.leftCols(numSearchDim);
        auto eigenvector = eigenvectors.col(eigenIndex);
        eigenvector = subSearchSpace * eigenSolver.getRawEigenvectors().col(0);
    }

    template<class ScalarType>
    size_t JacobiDavidson<ScalarType>::calcSearchSpaceDim(size_t order) {
        const size_t result = MaxBufferSize / order;
        if (result > MaxSearchDim)
            return MaxSearchDim;
        else if (result < MinSearchDim)
            return MinSearchDim;
        return result;
    }

    template<class ScalarType>
    inline void swap(JacobiDavidson<ScalarType>& __restrict obj1, JacobiDavidson<ScalarType>& __restrict obj2) noexcept {
        obj1.swap(obj2);
    }
}
