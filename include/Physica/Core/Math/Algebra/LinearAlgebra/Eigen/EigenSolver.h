/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/Schur.h"
#include "Physica/Core/Scalar/Complex.h" // IWYU pragma: export

namespace Physica {
    /**
     * A = XTX^{-1}
     * where X is matrix of eigenvectors
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     * [2] Eigen; https://eigen.tuxfamily.org/
     */
    template<Scalar T, size_t Order = Dynamic>
    class EigenSolver {
        using This = EigenSolver<T, Order>;
        constexpr static bool isComplex = T::isComplex;
    public:
        using RealType = T::RealType;
        using ComplexType = RealType::ComplexType;
        using EigenvalueVector = DenseVector<ComplexType, Order>;
        using EigenvectorMatrix = DenseMatrix<ComplexType, MatrixOption::Col | MatrixOption::Vector, Order, Order>;
        using RawEigenvectorType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
    private:
        using WorkingMatrix = Schur<T, Order>::WorkingMatrix;

        EigenvalueVector eigenvalues;
        RawEigenvectorType rawEigenvectors;
        bool computeEigenvectors;
    public:
        EigenSolver();
        EigenSolver(size_t size);
        template<Matrix M>
        EigenSolver(const M& source, bool computeEigenvectors_);
        EigenSolver(EigenvalueVector eigenvalues_, RawEigenvectorType rawEigenvectors_);
        EigenSolver(const This&) = default;
        EigenSolver(This&& solver) noexcept = default;
        ~EigenSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        void compute(const M& source, bool computeEigenvectors_);
        template<Matrix M>
        void compute_base(const M& source, bool computeEigenvectors_);
        template<Matrix M>
        void compute_mkl(const M& source, bool computeEigenvectors_);

        void sort();
        template<class Compare>
        void sort(Compare comp);
        void resize(size_t size);
        [[nodiscard]] auto reconstruct() const;
        [[nodiscard]] auto reconstruct_hermite() const;
        void swap(This& __restrict solver) noexcept;
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] EigenvalueVector& getEigenvalues() noexcept { return eigenvalues; }
        [[nodiscard]] EigenvectorMatrix getEigenvectors() const;
        /**
         * It is faster if all eigenvalues are real.
         */
        [[nodiscard]] const RawEigenvectorType& getRawEigenvectors() const noexcept { return rawEigenvectors; }
        /* Static members */
        static void sort(EigenvalueVector& eigenvalues, RawEigenvectorType& rawEigenvectors);
        template<class Compare>
        static void sort(EigenvalueVector& eigenvalues, RawEigenvectorType& rawEigenvectors, Compare comp);
    private:
        template<Matrix M>
        void pre_compute(const M& source, bool computeEigenvectors_);
        void computeRealMatEigenvalues(const WorkingMatrix& matrixT);
        void computeRealMatEigenvectors(WorkingMatrix& matrixT);
        void computeComplexMatEigenvectors(WorkingMatrix& matrixT);
        /* Static members */
        static bool defaultComp(ComplexType a, ComplexType b) noexcept;
    };

    template<Scalar T, size_t Order>
    EigenSolver<T, Order>::EigenSolver() : eigenvalues(), rawEigenvectors(), computeEigenvectors(false) {}

    template<Scalar T, size_t Order>
    EigenSolver<T, Order>::EigenSolver(size_t size) : EigenSolver() {
        resize(size);
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    EigenSolver<T, Order>::EigenSolver(const M& source, bool computeEigenvectors_)
            : computeEigenvectors(computeEigenvectors_) {
        resize(source.getRow());
        compute(source, computeEigenvectors);
    }
    /**
     * Reconstruct a matrix from its eigen decomposition
     */
    template<Scalar T, size_t Order>
    EigenSolver<T, Order>::EigenSolver(EigenvalueVector eigenvalues_, RawEigenvectorType rawEigenvectors_)
            : eigenvalues(std::move(eigenvalues_))
            , rawEigenvectors(std::move(rawEigenvectors_))
            , computeEigenvectors(true) {}

    template<Scalar T, size_t Order>
    template<Matrix M>
    void EigenSolver<T, Order>::compute(const M& source, bool computeEigenvectors_) {
        if constexpr (HasMKL())
            compute_mkl(source, computeEigenvectors_);
        else
            compute_base(source, computeEigenvectors_);
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    void EigenSolver<T, Order>::compute_base(const M& source, bool computeEigenvectors_) {
        pre_compute(source, computeEigenvectors_);
        if (source.getRow() == 1) [[unlikely]] {
            eigenvalues[0] = source.calc(0, 0);
            rawEigenvectors(0, 0) = T(1);
            return;
        }

        auto schur = Schur<T, Order>(source, computeEigenvectors);
        WorkingMatrix& matrixT = schur.getMatrixT();
        const size_t order = source.getRow();
        if constexpr (isComplex)
            eigenvalues = matrixT.diag();
        else
            computeRealMatEigenvalues(matrixT);

        if (computeEigenvectors) {
            if constexpr (isComplex)
                computeComplexMatEigenvectors(matrixT);
            else
                computeRealMatEigenvectors(matrixT);

            const auto& matrixU = schur.getMatrixU();
            for (size_t i = 0; i < order; ++i) {
                auto col = rawEigenvectors.col(i);
                col = matrixU.leftCols(i + 1) * matrixT.col(i).head(i + 1);
                if constexpr (isComplex)
                    col.toUnit();
            }
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::sort() {
        sort(defaultComp);
    }

    template<Scalar T, size_t Order>
    template<class Compare>
    void EigenSolver<T, Order>::sort(Compare comp) {
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
                rawEigenvectors.swap_col(i, index_min);
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::resize(size_t size) {
        assert((Order == Dynamic || Order == size) && "[Error]: size is not consistent");
        eigenvalues.resize(size);
        rawEigenvectors.resize(size, size);
    }

    template<Scalar T, size_t Order>
    auto EigenSolver<T, Order>::reconstruct() const {
        const size_t size = getSize();
        EigenvectorMatrix result(size, size);
        if constexpr (isComplex) {
            RawEigenvectorType inv = rawEigenvectors.inverse();
            for (size_t i = 0; i < size; ++i)
                result.col(i) = rawEigenvectors * hadamard(eigenvalues, inv.col(i));
            return result;
        }
        else {
            const EigenvectorMatrix eigenvectors = getEigenvectors();
            const EigenvectorMatrix inv = eigenvectors.inverse();
            for (size_t i = 0; i < size; ++i)
                result.col(i) = eigenvectors * hadamard(eigenvalues, inv.col(i));
            return RawEigenvectorType(result.reals());
        }
    }
    /**
     * Faster because no inverse matrix
     */
    template<Scalar T, size_t Order>
    auto EigenSolver<T, Order>::reconstruct_hermite() const {
        const size_t size = getSize();
        EigenvectorMatrix result(size, size);
        if constexpr (isComplex) {
            for (size_t i = 0; i < size; ++i)
                result.col(i) = rawEigenvectors * hadamard(eigenvalues, rawEigenvectors.conjugate().row(i));
            return result;
        }
        else {
            const EigenvectorMatrix eigenvectors = getEigenvectors();
            for (size_t i = 0; i < size; ++i)
                result.col(i) = eigenvectors * hadamard(eigenvalues, eigenvectors.conjugate().row(i));
            return RawEigenvectorType(result.reals());
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::swap(This& __restrict solver) noexcept {
        assert(this != &solver && "[Error]: Self swap is likely a bug");
        eigenvalues.swap(solver.eigenvalues);
        rawEigenvectors.swap(solver.rawEigenvectors);
        std::swap(computeEigenvectors, solver.computeEigenvectors);
    }

    template<Scalar T, size_t Order>
    EigenSolver<T, Order>::EigenvectorMatrix EigenSolver<T, Order>::getEigenvectors() const {
        assert(computeEigenvectors && "[Error]: Eigenvectors are not ready");
        if constexpr (isComplex)
            return rawEigenvectors;
        else {
            const size_t order = eigenvalues.getLength();
            EigenvectorMatrix result = EigenvectorMatrix(order, order);
            for (size_t i = 0; i < order; ++i) {
                if (eigenvalues[i].imag().isZero()) {
                    auto toCol = result.col(i);
                    toCol = rawEigenvectors.col(i);
                    toCol.toUnit();
                }
                else {
                    auto toCol1 = result.col(i);
                    auto toCol2 = result.col(i + 1);
                    auto fromCol1 = rawEigenvectors.col(i);
                    auto fromCol2 = rawEigenvectors.col(i + 1);
                    for (size_t j = 0; j < order; ++j) {
                        toCol1[j] = ComplexType(fromCol1[j].real(), fromCol2[j].real());
                        toCol2[j] = ComplexType(fromCol1[j].real(), -fromCol2[j].real());
                    }
                    toCol1.toUnit();
                    toCol2.toUnit();
                    ++i;
                }
            }
            return result;
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::sort(EigenvalueVector& eigenvalues, RawEigenvectorType& rawEigenvectors) {
        sort(eigenvalues, rawEigenvectors, defaultComp);
    }

    template<Scalar T, size_t Order>
    template<class Compare>
    void EigenSolver<T, Order>::sort(EigenvalueVector& eigenvalues, RawEigenvectorType& rawEigenvectors, Compare comp) {
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
            rawEigenvectors.swap_col(i, index_min);
        }
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    void EigenSolver<T, Order>::pre_compute([[maybe_unused]] const M& source, bool computeEigenvectors_) {
        static_assert(std::is_same<T, typename M::ScalarType>::value, "[Error]: Inconsistent ScalarType");
        assert(source.getRow() > 0);
        assert(source.getRow() == source.getCol());
        assert(source.getRow() == eigenvalues.getLength());
        computeEigenvectors = computeEigenvectors_;
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::computeRealMatEigenvalues(const WorkingMatrix& matrixT) {
        const size_t order = matrixT.getRow();
        for (size_t i = 0; i < order;) {
            if (i == order - 1 || matrixT(i + 1, i).isZero()) {
                eigenvalues[i] = matrixT(i, i);
                i += 1;
            }
            else {
                const T p = T(0.5) * (matrixT(i, i) - matrixT(i + 1, i + 1));
                T z;
                /* Referenced from eigen, to avoid overflow */ {
                    T t0 = matrixT(i + 1, i);
                    T t1 = matrixT(i, i + 1);
                    const T max = std::max(abs(p), std::max(abs(t0), abs(t1)));
                    const T inv_max = reciprocal(max);
                    t0 *= inv_max;
                    t1 *= inv_max;
                    T p0 = p * inv_max;
                    z = max * sqrt(abs(square(p0) + t0 * t1));
                }
                const T real = p + matrixT(i + 1, i + 1);
                eigenvalues[i] = Complex<T>(real, z);
                eigenvalues[i + 1] = Complex<T>(real, -z);
                i += 2;
            }
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::computeRealMatEigenvectors(WorkingMatrix& matrixT) {
        const size_t order = matrixT.getRow();
        for (size_t i = order - 1; i < order; --i) {
            auto block = matrixT.topLeftCorner(i + 1, i + 1);
            if (eigenvalues[i].imag().isZero()) {
                auto col = block.col(i);
                col[i] = T(1);
                for (size_t j = i - 1; j < i; --j) {
                    const T deltaEigen = eigenvalues[i].real() - eigenvalues[j].real();
                    if (eigenvalues[j].imag().isZero()) {
                        const T factor = deltaEigen.isZero() ? T(std::numeric_limits<T>::epsilon()) : deltaEigen;
                        col[j] = (col.tail(j + 1) * block.row(j).tail(j + 1)) / factor;
                    }
                    else {
                        const auto row1 = block.row(j - 1);
                        const auto row2 = block.row(j);
                        const auto tail = col.tail(j + 1);
                        const T dot1 = tail * row1.tail(j + 1);
                        const T dot2 = tail * row2.tail(j + 1);
                        const T inv_determinate = reciprocal(square(deltaEigen) + square(eigenvalues[j].imag()));
                        col[j - 1] = (dot1 * block(j, j) - dot2 * block(j - 1, j)) * inv_determinate;
                        col[j] = (dot2 * block(j - 1, j - 1) - dot1 * block(j, j - 1)) * inv_determinate;
                        --j;
                    }
                }
            }
            else {
                assert(eigenvalues[i].imag().isNegative());
                auto col1 = block.col(i - 1);
                auto col2 = block.col(i);
                // Referenced from eigen, to ensure numerical stable
                if (abs(matrixT(i, i - 1)) > abs(matrixT(i - 1, i))) {
                    auto temp = reciprocal(matrixT(i, i - 1));
                    col1[i - 1] = eigenvalues[i].imag() * temp;
                    col2[i - 1] = (eigenvalues[i].real() - matrixT(i, i)) * temp;
                }
                else {
                    const auto c = Complex<T>(T(0), -matrixT(i - 1, i)) / Complex<T>(matrixT(i - 1, i - 1) - eigenvalues[i].real(), eigenvalues[i].imag());
                    col1[i - 1] = c.real();
                    col2[i - 1] = c.imag();
                }
                col1[i] = T(0);
                col2[i] = T(1);

                for (size_t j = i - 2; j < i; --j) {
                    if (eigenvalues[j].imag().isZero()) {
                        auto row = block.row(j);
                        auto tail = row.tail(j + 1);
                        const T dot1 = -(tail * col1.tail(j + 1));
                        const T dot2 = -(tail * col2.tail(j + 1));
                        const T a = block(j, j) - eigenvalues[i].real();
                        const T b = eigenvalues[i].imag();
                        const T inv_denominator = reciprocal(square(a) + square(b));
                        col1[j] = (a * dot1 + b * dot2) * inv_denominator;
                        col2[j] = (a * dot2 - b * dot1) * inv_denominator;
                    }
                    else {
                        auto row1 = block.row(j - 1);
                        auto tail1 = row1.tail(j + 1);
                        const T dot11 = tail1 * col1.tail(j + 1);
                        const T dot12 = tail1 * col2.tail(j + 1);
                        auto row2 = block.row(j);
                        auto tail2 = row2.tail(j + 1);
                        const T dot21 = tail2 * col1.tail(j + 1);
                        const T dot22 = tail2 * col2.tail(j + 1);

                        const T x = matrixT(j - 1, j);
                        const T y = matrixT(j, j - 1);
                        const T deltaEigen = eigenvalues[j].real() - eigenvalues[i].real();
                        const T vr = square(deltaEigen) + square(eigenvalues[j].imag()) - square(eigenvalues[i].imag());
                        const T vi = T(2) * deltaEigen * eigenvalues[i].imag();

                        Complex<T> c{};
                        const T temp1 = matrixT(j, j) - eigenvalues[i].real();
                        c.real() = x * dot21 - temp1 * dot11 + eigenvalues[i].imag() * dot12;
                        c.imag() = x * dot22 - temp1 * dot12 - eigenvalues[i].imag() * dot11;
                        c /= Complex<T>(vr, vi);
                        matrixT(j - 1, i - 1) = c.real();
                        matrixT(j - 1, i) = c.imag();

                        if (abs(x) > (abs(temp1) + abs(eigenvalues[i].imag()))) {
                            const T temp2 = matrixT(j - 1, j - 1) - eigenvalues[i].real();
                            matrixT(j, i - 1) = (-dot11 - temp2 * matrixT(j - 1, i - 1) + eigenvalues[i].imag() * matrixT(j - 1, i)) / x;
                            matrixT(j, i) = (-dot12 - temp2 * matrixT(j - 1, i) - eigenvalues[i].imag() * matrixT(j - 1, i - 1)) / x;
                        }
                        else {
                            c = Complex<T>(-dot21 - y * matrixT(j - 1, i - 1), -dot22 - y * matrixT(j - 1, i)) / Complex<T>(temp1, eigenvalues[i].imag());
                            matrixT(j, i - 1) = c.real();
                            matrixT(j, i) = c.imag();
                        }
                        --j;
                    }
                }
                --i;
            }
        }
    }

    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::computeComplexMatEigenvectors(WorkingMatrix& matrixT) {
        const size_t order = matrixT.getRow();
        for (size_t i = order - 1; i < order; --i) {
            auto block = matrixT.topLeftCorner(i + 1, i + 1);
            auto col = block.col(i);
            col[i] = T(1);
            for (size_t j = i - 1; j < i; --j) {
                assert(!scalarNear(eigenvalues[i], eigenvalues[j], 1E-14)); // TODO: handle degeneracy
                auto row = block.row(j);
                col[j] = (col.tail(j + 1) * row.tail(j + 1)) / (eigenvalues[i] - eigenvalues[j]);
            }
        }
    }

    template<Scalar T, size_t Order>
    bool EigenSolver<T, Order>::defaultComp(ComplexType a, ComplexType b) noexcept {
        return a.real() < b.real();
    }
}

namespace std {
    template<Physica::Scalar T, size_t Order>
    inline void swap(Physica::EigenSolver<T, Order>& __restrict solver1,
                     Physica::EigenSolver<T, Order>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}

#ifdef PHYSICA_MKL
    #include "EigenSolver_MKL.h" // IWYU pragma: export
#endif
