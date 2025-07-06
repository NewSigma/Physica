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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    /**
     * Decompose matrix like A = QR(no pivoting), or A = QRP(poviting)
     */
    template<Scalar T, bool Pivot = false>
    class QRDecomp {
        constexpr static bool isComplex = T::isComplex;
        using This = QRDecomp<T, Pivot>;
        using MatrixND = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
        using Tc = T::ComplexType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
        using Perm = std::conditional<Pivot, PermMatrix<T>, PlainStruct<void>>::type;
    private:
        MatrixND working;
        VectorND<T> taus;
        [[no_unique_address]] Perm perm;
    public:
        QRDecomp() = default;
        QRDecomp(size_t row, size_t col);
        QRDecomp(const Matrix auto& source);
        QRDecomp(const This&) = default;
        QRDecomp(This&&) noexcept = default;
        ~QRDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        void compute(const M& source);
        void compute_base(const Matrix auto& source);
        void compute_mkl(const Matrix auto& source);

        [[nodiscard]] T calcDetQ() const;
        void toQDT(VectorND<T>& diagD);
        [[nodiscard]] VectorND<T> toQDT();

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto& getWorking() noexcept { return working; }
        [[nodiscard]] const auto& getWorking() const noexcept { return working; }
        [[nodiscard]] const auto& getTaus() const noexcept { return taus; }
        [[nodiscard]] size_t getRow() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return working.getCol(); }
        [[nodiscard]] MatrixND getMatrixQ() const;
        [[nodiscard]] auto getMatrixR() const noexcept;
        [[nodiscard]] MatrixND getMatrixQ_mkl() const;
        [[nodiscard]] const auto& getMatrixP() const noexcept requires(Pivot) { return perm; }
    };

    template<Scalar T, bool Pivot>
    QRDecomp<T, Pivot>::QRDecomp(size_t row, size_t col) {
        resize(row, col);
    }

    template<Scalar T, bool Pivot>
    QRDecomp<T, Pivot>::QRDecomp(const Matrix auto& source) : QRDecomp(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Scalar T, bool Pivot>
    template<Matrix M>
    void QRDecomp<T, Pivot>::compute(const M& source) {
        assert(getRow() == source.getRow());
        assert(getCol() == source.getCol());
        if constexpr (HasMKL() && (M::SizeAtCompile > 16 || M::SizeAtCompile == Dynamic))
            compute_mkl<M>(source);
        else
            compute_base<M>(source);
    }

    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::compute_base(const Matrix auto& source) {
        assert(!Pivot && "Not implemented");
        working = source;

        size_t i = 0;
        for (; i < taus.getLength() - 1; ++i) {
            auto col = working.col(i);
            auto buffer = col.tail(i);
            const auto unit = buffer[0].unit();
            const auto norm = buffer.householder();
            auto corner = working.bottomRightCorner(i, i + 1);
            applyHouseholder(buffer, corner);
            taus[i] = std::exchange(col[i], -norm * unit);
        }
    }

    template<Scalar T, bool Pivot>
    T QRDecomp<T, Pivot>::calcDetQ() const {
        int sign = 0;
        for (auto tau : taus)
            sign += tau.isPositive();
        return T(sign % 2 == 0 ? 1.0 : -1.0);
    }
    /**
     * Decompose matrix like A = QDT(no pivoting), or A = QDTP(poviting), where D is diagonal
     */
    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::toQDT(VectorND<T>& diagD) {
        const size_t length = taus.getLength();
        assert(diagD.getLength() == length);
        for (size_t i = 0; i < length; ++i) {
            if (working(i, i).isZero()) {
                diagD[i] = 1;
                continue;
            }
            diagD[i] = working(i, i);
            working.row(i).tail(i) *= reciprocal(diagD[i]);
        }
    }

    template<Scalar T, bool Pivot>
    VectorND<T> QRDecomp<T, Pivot>::toQDT() {
        VectorND<T> vecD(taus.getLength());
        toQDT(vecD);
        return vecD;
    }

    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::resize(size_t row, size_t col) {
        working.resize(row, col);
        auto l = std::min(row, col);
        taus.resize(l);
        taus[l - 1] = 0; // For historic reason, BLAS-like interface will allocate a unused element
        if constexpr (Pivot)
            perm = Perm(col);
    }

    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        taus.swap(obj.taus);
        perm.swap(obj.perm);
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::getMatrixQ() const -> MatrixND {
        if constexpr (HasMKL())
            return getMatrixQ_mkl();
        else
            noImpl(__func__);
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::getMatrixR() const noexcept {
        return working.triu();
    }
}

#ifdef PHYSICA_MKL
    #include "QRDecomp_MKL.h"
#endif
