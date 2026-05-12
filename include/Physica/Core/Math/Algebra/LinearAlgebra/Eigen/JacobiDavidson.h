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

#include "EigenSolver.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearSystem/LinearCG.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Orthogonalize.h"

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
        using DotProjectSpace = std::conditional<T::isComplex(), DenseHermiteMatrix<T>, DenseSymmMatrix<T>>::type;

        constexpr static size_t MaxIterationPerEigen = 4;
        constexpr static size_t MaxSearchDim = 32; // Refer to [1]
        constexpr static size_t MinSearchDim = 5; // Refer to [1]
        constexpr static size_t MaxBufferSize = 1024UL * 1024UL * 1024UL / sizeof(T);
        constexpr static size_t MaxLinearSolverIteration = 64;
        constexpr static double LinearSolverPrecision = 1E-4;
        constexpr static double DefaultStableThreshold = 1E-5;
        constexpr static double NearConvergeThreshold = 1E-3;
    public:
        constexpr static Tr InvalidGoal = std::numeric_limits<Tr>::max();
    private:
        LinearCG<T> linearSolver;
        EigenSolver<T> eigenSolver;

        size_t curSearchDim;
        MatrixND<T> searchSpace;
        MatrixND<T> dotSpace;
        MatrixND<T> searchSpaceProj;
        DotProjectSpace dProj;
        MatrixND<T> eigenBuffer;

        VectorND<T> eigenvalues;
        MatrixND<T> eigenvectors;

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
        void compute(const Matrix auto& source, const VectorND<T>& initial, T eigenGoal = T(InvalidGoal));
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
        void initSearchSpace(const Matrix auto& source, const VectorND<T>& initial);
        void projSearchSpace(const Matrix auto& source, size_t index);
        void expandSpace(const Matrix auto& source, size_t dim);

        [[nodiscard]] bool findEigenPair(const Matrix auto& source, size_t index, VectorND<T>& residual, VectorND<T>& buffer, T& eigenGoal);
        [[nodiscard]] T updateGoal(const Matrix auto& source, VectorND<T>& buffer, T eigenvalue, Tr squaredRes, Tr& lastDeltaEigen, T eigenGoal);
        void extremeSearch();
        void refinedSearch(size_t index, T eigenGoal);
        void correction(size_t index, T eigenGoal, const Matrix auto& source, const VectorND<T>& residual, VectorND<T>& buffer);
        /* Static member */
        static size_t calcSearchSpaceDim(size_t order) noexcept;
    };

    template<Scalar T>
    JacobiDavidson<T>::JacobiDavidson() : curSearchDim() {
        linearSolver.mustConverge = false;
        linearSolver.setError(LinearSolverPrecision);
        linearSolver.setIterationLimit(MaxLinearSolverIteration);
    }

    template<Scalar T>
    JacobiDavidson<T>::JacobiDavidson(size_t size, size_t numRequired) : JacobiDavidson() {
        linearSolver.resize(size);
        resize(size, numRequired);
    }

    template<Scalar T>
    void JacobiDavidson<T>::compute(const Matrix auto& source, const VectorND<T>& initial, T eigenGoal) {
        static_assert(source.isStaticHermite(), "[Error]: Support for complex eigenvalues is not implemented");
        assert(source.getRow() == source.getCol() && "[Error]: Matrix should be square");
        assert(source.getRow() == initial.getLength() && "[Error]: Dimensions do not match");
        assert(eigenGoal == InvalidGoal && "[Error]: Not implemented");
        initSearchSpace(source, initial);

        VectorND<T> residual(initial.getLength());
        VectorND<T> buffer(initial.getLength());
        for (size_t i = 0; i < getNumRequired(); ++i) {
            if (i > 0)
                projSearchSpace(source, i);

            size_t iteration = 0;
            while (!findEigenPair(source, i, residual, buffer, eigenGoal)) {
                iteration += 1;
                if (iteration == MaxIterationPerEigen) [[unlikely]]
                    throw BadConvergenceException("Exceed max iteration of JacobiDavidson");

                auto v = eigenvectors.col(i);
                residual = source * v;
                eigenvalues[i] = v.conjugate() * residual;
                eigenGoal = eigenvalues[i];
                initSearchSpace(source, v);
                for (size_t j = 1; j <= i; ++j)
                    projSearchSpace(source, j);
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
        assert(size > MaxSearchDim && "[Error]: Please use dense eigen solver for small matrix");
        const size_t dim = calcSearchSpaceDim(size);
        linearSolver.resize(size);
        eigenvalues.resize(numRequired);
        eigenvectors.resize(size, numRequired);
        eigenSolver.resize(dim, true);
        searchSpace.resize(size, dim);
        dotSpace.resize(size, dim);
        searchSpaceProj.resize(dim, dim);
        dProj.resize(dim, dim);
        eigenBuffer.resize(dim, dim);
    }

    template<Scalar T>
    void JacobiDavidson<T>::swap(JacobiDavidson<T>& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        linearSolver.swap(obj.linearSolver);
        eigenSolver.swap(obj.eigenSolver);

        std::swap(curSearchDim, obj.curSearchDim);
        searchSpace.swap(obj.searchSpace);
        dotSpace.swap(obj.dotSpace);
        searchSpaceProj.swap(obj.searchSpaceProj);
        dProj.swap(obj.dProj);
        eigenBuffer.swap(obj.eigenBuffer);

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
    void JacobiDavidson<T>::initSearchSpace(const Matrix auto& source, const VectorND<T>& initial) {
        /* dim = 0 */ {
            auto v = searchSpace.col(0);
            v = initial;
            v.toUnit();
            expandSpace(source, 0);

            auto buffer = searchSpace.col(1);
            buffer = dotSpace.col(0) - searchSpaceProj[0, 0] * v;
            const Tr squaredRes = buffer.squaredNorm();
            const bool isGoodInitial = squaredRes < Tr(NearConvergeThreshold);
            if (isGoodInitial)
                curSearchDim = 1;
        }
        for (size_t dim = 1; dim < MinSearchDim; ++dim) {
            searchSpace.col(dim) = source * searchSpace.col(dim - 1);
            expandSpace(source, dim);
        }
        curSearchDim = MinSearchDim;
    }
    /**
     * New eigen-pair found; project it out of the current searching space and restart
     * 
     * Refer to [1]
     */
    template<Scalar T>
    void JacobiDavidson<T>::projSearchSpace(const Matrix auto& source, size_t index) {
        assert(index > 0);
        const auto lastEigenvector = eigenvectors.col(index - 1);
        for (size_t dim = 0; dim < MinSearchDim; ++dim) {
            auto v = searchSpace.col(dim);
            T dot = lastEigenvector.conjugate() * v;
            Tr normalizer = reciprocal(sqrt(Tr(1) - dot.squaredNorm()));
            v = (v - dot * lastEigenvector) * normalizer;
            expandSpace(source, dim);
        }
        curSearchDim = MinSearchDim;
    }

    template<Scalar T>
    void JacobiDavidson<T>::expandSpace(const Matrix auto& source, size_t dim) {
        const size_t i = dim;
        auto dir = searchSpace.col(i);
        if (i != 0)
            gram_schmidt(searchSpace.leftCols(i), dir, searchSpaceProj.col(i).head(i));

        auto dot = dotSpace.col(i);
        dot = source * dir;
        if (i != 0) {
            auto v = searchSpaceProj.col(i).head(i);
            v = dotSpace.leftCols(i).hermite() * dir;
            if constexpr (source.isStaticHermite())
                searchSpaceProj.row(i).head(i) = v.conjugate();
            else
                noImpl();

            for (size_t j = 0; j < i; ++j)
                dProj[j, i] = dotSpace.col(j).conjugate() * dot;
        }
        searchSpaceProj[i, i] = dot.conjugate() * dir;
        dProj[i, i] = dot.squaredNorm();
    }

    template<Scalar T>
    bool JacobiDavidson<T>::findEigenPair(const Matrix auto& source, size_t index, VectorND<T>& residual, VectorND<T>& buffer, T& eigenGoal) {
        Tr lastDeltaEigen = std::numeric_limits<Tr>::max();
        T& eigenvalue = eigenvalues[index];
        auto eigenvector = eigenvectors.col(index);
        while (true) {
            const bool isGoodInitial = curSearchDim == 1;
            if (isGoodInitial) {
                eigenvector = searchSpace.col(0);
                residual = searchSpace.col(1);
                eigenvalue = searchSpaceProj[0, 0];
                eigenGoal = eigenvalue;
            }
            else {
                bool isFirstEigen = index == 0;
                if (isFirstEigen) {
                    extremeSearch();
                    eigenvalue = eigenSolver.getEigenvalues()[0].real();
                    residual = source * eigenvector - eigenvalue * eigenvector;
                    eigenGoal = std::min(eigenGoal.real(), eigenvalue.real());
                }
                else {
                    refinedSearch(index, eigenGoal);
                    residual = source * eigenvector;
                    eigenvalue = eigenvector.conjugate() * residual;
                    residual -= eigenvalue * eigenvector;
                }
                Tr squaredRes = residual.squaredNorm();
                bool converged = squaredRes < error;
                bool shouldRestart = curSearchDim == searchSpace.getCol();
                if (converged || shouldRestart)
                    return converged;

                if (!isFirstEigen)
                    eigenGoal = updateGoal(source, buffer, eigenvalue, squaredRes, lastDeltaEigen, eigenGoal);
            }
            correction(index, eigenGoal, source, residual, buffer);
        }
    }

    template<Scalar T>
    T JacobiDavidson<T>::updateGoal(const Matrix auto& source, VectorND<T>& buffer, T eigenvalue, Tr squaredRes, Tr& lastDeltaEigen, const T eigenGoal) {
        const VectorND<T> eigenvector2 = searchSpace.leftCols(curSearchDim) * eigenSolver.getRawEigenvectors().col(1);
        buffer = source * eigenvector2;
        const T eigenvalue2 = eigenvector2.conjugate() * buffer;
        const Tr deltaEigen = abs(eigenvalue2.real() - eigenvalue.real());
        const bool nearConverge = squaredRes < square(deltaEigen);
        const bool deltaEigenStable = abs(deltaEigen - lastDeltaEigen) < stableThreshold * lastDeltaEigen;
        const bool goodDeltaEigen = lastDeltaEigen > std::numeric_limits<Tr>::epsilon();
        const bool increaseGoal = !goodDeltaEigen || (deltaEigenStable && nearConverge);
        const bool decreaseGoal = eigenvalue.real() < eigenGoal.real();
        if (increaseGoal || decreaseGoal)
            return eigenvalue;
        lastDeltaEigen = deltaEigen;
        return eigenGoal;
    }

    template<Scalar T>
    void JacobiDavidson<T>::extremeSearch() {
        assert(curSearchDim > 1 && "[Error]: No need to search if dim = 1");
        eigenSolver.resize(curSearchDim);
        eigenSolver.compute(searchSpaceProj.topLeftCorner(curSearchDim));
        eigenSolver.sort();
        eigenvectors.col(0) = searchSpace.leftCols(curSearchDim) * eigenSolver.getRawEigenvectors().col(0);
    }

    template<Scalar T>
    void JacobiDavidson<T>::refinedSearch(size_t index, T eigenGoal) {
        assert(curSearchDim > 1 && "[Error]: No need to search if dim = 1");
        assert(eigenGoal.real() != Tr(InvalidGoal) && "[Error]: Goal is not updated, this is a bug");
        auto corner = eigenBuffer.topLeftCorner(curSearchDim);
        {
            const auto proj = dProj.topLeftCorner(curSearchDim);
            const auto basis = searchSpaceProj.topLeftCorner(curSearchDim);
            if constexpr (T::isComplex())
                corner = proj - eigenGoal * basis - eigenGoal.conjugate() * basis.hermite();
            else
                corner = proj - T(2) * eigenGoal * basis;
        }
        eigenSolver.resize(curSearchDim);
        eigenSolver.compute(corner);
        eigenSolver.sort();
        eigenvectors.col(index) = searchSpace.leftCols(curSearchDim) * eigenSolver.getRawEigenvectors().col(0);
    }

    template<Scalar T>
    void JacobiDavidson<T>::correction(size_t index, T eigenGoal, const Matrix auto& source, const VectorND<T>& residual, VectorND<T>& buffer) {
        const auto orthogonalSpace = eigenvectors.leftCols(index + 1);
        auto corrector = searchSpace.col(curSearchDim);
        corrector = residual;
        gram_schmidt(orthogonalSpace, corrector);
        linearSolver.cg([this, index, eigenGoal, &buffer, &source](const VectorND<T>& v, VectorND<T>& dot) {
            const size_t dim = index + 1;
            const auto orthogonalSpace = eigenvectors.leftCols(dim);
            {
                auto head = dot.head(dim);
                head = orthogonalSpace.hermite() * v;
                buffer = v - orthogonalSpace * head;
            }
            dot = source * buffer - eigenGoal * buffer;
            {
                auto head = buffer.head(dim);
                head = orthogonalSpace.hermite() * dot;
                dot -= orthogonalSpace * head;
            }
        }, corrector);
        expandSpace(source, curSearchDim);
        curSearchDim += 1;
    }

    template<Scalar T>
    size_t JacobiDavidson<T>::calcSearchSpaceDim(size_t order) noexcept {
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
