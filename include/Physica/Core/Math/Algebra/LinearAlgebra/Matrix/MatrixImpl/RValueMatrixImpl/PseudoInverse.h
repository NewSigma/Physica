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

#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    class SVD;

    template<Scalar T, size_t Order>
    class DiagMatrix;

    template<Matrix M>
    class PseudoInverse : public RValueMatrix<PseudoInverse<M>> {
        using This = PseudoInverse<M>;
        using Base = RValueMatrix<This>;

        using typename Base::T;

        const M& mat;
    public:
        explicit PseudoInverse(const M& mat_);
        PseudoInverse(const This&) = default;
        PseudoInverse(This&&) noexcept = default;
        ~PseudoInverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const M& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
        [[nodiscard]] size_t getOrder() const noexcept { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M>
    PseudoInverse<M>::PseudoInverse(const M& mat_) : mat(mat_) {}

    template<Matrix M>
    void PseudoInverse<M>::assign(Matrix auto& target) const {
        SVD<T, mat.getRowAtCompile(), mat.getColAtCompile()> svd(mat);
        auto diagD = DiagMatrix<T, decltype(svd)::NumSingularValue>(svd.getNumSingular());
        const auto& singular = svd.getSingulars();
        const T tol = singular.max() * std::numeric_limits<T>::epsilon();
        for (size_t i = 0; i < singular.getLength(); ++i)
            diagD.diag()[i] = (singular.calc(i) < tol) ? T(0) : reciprocal(singular.calc(i));
        target = svd.getMatrixV() * diagD * svd.getMatrixU().hermite();
    }

    template<Matrix M>
    __host__ __device__ consteval bool PseudoInverse<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<PseudoInverse<M>> : public Traits<Transpose<M>> {};
}
