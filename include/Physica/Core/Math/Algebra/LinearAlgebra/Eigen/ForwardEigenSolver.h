/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "EigenSolver.h"

namespace Physica {
    template<Scalar Tv, size_t Order>
    class EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order> {
        static_assert(!Tv::isDiffable, "[Error]: Invalid tparam");
        using T = Diff<Tv, DiffMode::Forward, 1>;
        using This = EigenSolver<T, Order>;

        using Tr = T::RealType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        using EigenvalueVector = DenseVector<Tc, Order>;
        using EigenvectorMatrix = DenseMatrix<Tc, MatrixMajor::Col, Order, Order>;
    private:
        EigenvalueVector eigenvalues;
        EigenvectorMatrix eigenvectors;
    public:
        EigenSolver() = default;
        explicit EigenSolver(size_t size, bool needEigenvectors);
        explicit EigenSolver(const Matrix auto& source, bool needEigenvectors);
        EigenSolver(const This&) = default;
        EigenSolver(This&&) noexcept = default;
        ~EigenSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);

        void sort();
        void sort(std::invocable<Tcv, Tcv> auto comp);
        void resize(size_t size, bool needEigenvectors);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return eigenvalues.getLength(); }
        [[nodiscard]] auto&& getEigenvalues(this auto&& self) noexcept;
        [[nodiscard]] const auto& getEigenvectors() const;
        [[nodiscard]] bool getNeedEigenvectors() const noexcept { return !eigenvectors.empty(); }
    };

    template<Scalar Tv, size_t Order>
    EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::EigenSolver(size_t size, bool needEigenvectors) {
        resize(size, needEigenvectors);
    }

    template<Scalar Tv, size_t Order>
    EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::EigenSolver(const Matrix auto& source, bool needEigenvectors)
            : EigenSolver(source.getRow(), needEigenvectors) {
        compute(source);
    }

    template<Scalar Tv, size_t Order>
    void EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::compute(const Matrix auto& source) {
        using GradMatrix = EigenSolver<Tv, Order>::EigenvectorMatrix;
        EigenSolver<Tv, Order> solver(source.values(), true);
        const auto basis = solver.getEigenvectors();

        GradMatrix transGrads;
        constexpr bool isHermite = source.isStaticHermite();
        constexpr bool isRealSymm = !T::isComplex && source.isStaticSymm();
        if constexpr (isHermite || isRealSymm)
            transGrads = basis.hermite() * GradMatrix(source.grads() * basis);
        else
            transGrads = GradMatrix(basis.inv()) * GradMatrix(source.grads() * basis);

        eigenvalues.values() = std::move(solver.getEigenvalues());
        eigenvalues.grads() = transGrads.diag();

        if (getNeedEigenvectors()) {
            const double threshold = std::sqrt(std::numeric_limits<T>::epsilon());
            size_t size = getSize();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < size; ++j) {
                    bool isDegenerate = scalarNear(eigenvalues.calc_value(i), eigenvalues.calc_value(j), threshold);
                    if (isDegenerate)
                        transGrads[i, j] = 0;
                    else
                        transGrads[i, j] /= (eigenvalues.calc_value(j) - eigenvalues.calc_value(i));
                }
            }
            eigenvectors.values() = std::move(basis);
            eigenvectors.grads() = basis * transGrads;
        }
    }

    template<Scalar Tv, size_t Order>
    void EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::sort() {
        sort(EigenSolver<Tv, Order>::defaultComp);
    }

    template<Scalar Tv, size_t Order>
    void EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::sort(std::invocable<Tcv, Tcv> auto comp) {
        const size_t order = eigenvalues.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (comp(eigenvalues[j].value(), eigenvalues[index_min].value()))
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

    template<Scalar Tv, size_t Order>
    void EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::resize(size_t size, bool needEigenvectors) {
        assert((Order == Dynamic || Order == size) && "[Error]: size is not consistent");
        eigenvalues.resize(size);
        if (needEigenvectors)
            eigenvectors.resize(size, size);
    }

    template<Scalar Tv, size_t Order>
    void EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        eigenvalues.swap(obj.eigenvalues);
        eigenvectors.swap(obj.eigenvectors);
    }

    template<Scalar Tv, size_t Order>
    auto&& EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::getEigenvalues(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).eigenvalues;
    }

    template<Scalar Tv, size_t Order>
    const auto& EigenSolver<Diff<Tv, DiffMode::Forward, 1>, Order>::getEigenvectors() const {
        assert(getNeedEigenvectors() && "[Error]: Eigenvectors are not ready");
        return eigenvectors;
    }
}
