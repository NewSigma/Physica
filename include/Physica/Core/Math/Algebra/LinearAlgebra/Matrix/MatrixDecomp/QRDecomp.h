/*
 * Copyright 2025 Weibo He.
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

namespace Physica {
    template<Scalar T>
    class QRDecomp {
        using This = QRDecomp<T>;
        using MatrixND = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
        using Tc = T::ComplexType;
        constexpr static bool isComplex = T::isComplex;
    private:
        MatrixND working;
        VectorND<T> taus;
    public:
        QRDecomp() = default;
        QRDecomp(size_t row, size_t col);
        template<Matrix M>
        QRDecomp(const M& source);
        QRDecomp(const This&) = default;
        QRDecomp(This&&) noexcept = default;
        ~QRDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        void compute(const M& source, bool pivote = false);
        template<Matrix M>
        void compute_mkl(const M& source, bool pivote = false);

        [[nodiscard]] VectorND<T> toQTD();

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getWorking() const noexcept { return working; }
        [[nodiscard]] const auto& getTaus() const noexcept { return taus; }
        [[nodiscard]] size_t getRow() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return working.getCol(); }
        [[nodiscard]] MatrixND getMatrixQ() const noexcept;
        [[nodiscard]] auto getMatrixR() const noexcept;
        [[nodiscard]] MatrixND getMatrixQ_mkl() const noexcept;
    };

    template<Scalar T>
    QRDecomp<T>::QRDecomp(size_t row, size_t col) {
        resize(row, col);
    }

    template<Scalar T>
    template<Matrix M>
    QRDecomp<T>::QRDecomp(const M& source) : QRDecomp(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Scalar T>
    template<Matrix M>
    void QRDecomp<T>::compute(const M& source, bool pivote) {
        assert(getRow() == source.getRow());
        assert(getCol() == source.getCol());
        if constexpr (HasMKL())
            compute_mkl(source, pivote);
        else
            noImpl(__func__);
    }

    template<Scalar T>
    VectorND<T> QRDecomp<T>::toQTD() {
        const size_t length = taus.getLength();
        VectorND<T> diag(length);
        for (size_t i = 0; i < length; ++i) {
            diag[i] = working(i, i);
            working.col(i).head(i + 1) *= reciprocal(diag[i]);
        }
        return diag;
    }

    template<Scalar T>
    void QRDecomp<T>::resize(size_t row, size_t col) {
        working.resize(row, col);
        taus.resize(std::min(row, col));
    }

    template<Scalar T>
    void QRDecomp<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj);
        taus.swap(obj.taus);
    }

    template<Scalar T>
    auto QRDecomp<T>::getMatrixQ() const noexcept -> MatrixND {
        if constexpr (HasMKL())
            return getMatrixQ_mkl();
        else
            noImpl(__func__);
    }

    template<Scalar T>
    auto QRDecomp<T>::getMatrixR() const noexcept {
        return working.triu();
    }
}

#ifdef PHYSICA_MKL
    #include "QRDecomp_MKL.h"
#endif
