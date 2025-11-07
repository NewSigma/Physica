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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    /**
     * Decompose matrix like A = QR(no pivoting), or A = QRP(poviting)
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:249
     */
    template<Scalar T, bool Pivot = false>
    class QRDecomp {
        constexpr static bool isComplex = T::isComplex;
        using This = QRDecomp<T, Pivot>;
        using Perm = std::conditional<Pivot, PermMatrix<T>, PlainStruct<void>>::type;

        using Tr = T::RealType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tv = T::ValueType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
        constexpr static size_t Threhold = 16;
    private:
        MatrixND<T> working;
        VectorND<T> taus;
        [[no_unique_address]] Perm perm;
    public:
        QRDecomp() = default;
        explicit QRDecomp(size_t order);
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

        [[nodiscard]] Tv calcDetQ() const noexcept;
        void toQDT(VectorND<Tr>& diagD) noexcept;
        [[nodiscard]] VectorND<Tr> toQDT();

        [[nodiscard]] T det() const noexcept;

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto& getWorking() noexcept { return working; }
        [[nodiscard]] const auto& getWorking() const noexcept { return working; }
        [[nodiscard]] const auto& getTaus() const noexcept { return taus; }
        [[nodiscard]] size_t getRow() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return working.getCol(); }
        [[nodiscard]] MatrixND<T> getMatrixQ() const;
        [[nodiscard]] MatrixND<T> getMatrixQ_mkl() const;
        [[nodiscard]] MatrixND<T> getMatrixQ_base() const;
        [[nodiscard]] auto getMatrixR() const noexcept;
        [[nodiscard]] const auto& getMatrixP() const noexcept requires(Pivot) { return perm; }
    };

    template<Scalar T, bool Pivot>
    QRDecomp<T, Pivot>::QRDecomp(size_t order) : QRDecomp(order , order) {}

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
        constexpr bool SmallMatrix = M::SizeAtCompile <= Threhold && M::SizeAtCompile != Dynamic;
        if constexpr (HasMKL() && !SmallMatrix) {
            if constexpr (M::SizeAtCompile > Threhold) {
                compute_mkl<M>(source);
                return;
            }

            if (source.getSize() > Threhold)
                compute_mkl<M>(source);
            else
                compute_base<M>(source);
        }
        else
            compute_base<M>(source);
    }

    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::compute_base(const Matrix auto& source) {
        static_assert(!Pivot && "Not implemented");
        working = source;

        size_t i = 0;
        for (; i < taus.getLength() - !working.isOverdetermined(); ++i) {
            auto col = working.col(i);
            auto buffer = col.tail(i);
            const T sign = unit(buffer[0]);
            const Tr norm = buffer.householder();
            bool isFinalColumn = i + 1 >= getCol();
            if (!isFinalColumn)
                applyHouseholder(buffer, working.bottomRightCorner(i, i + 1));
            taus[i] = std::exchange(col[i], -norm * sign);
        }
        // Other LAPACK implementations might modify the final element, so always clear it.
        if (i < taus.getLength())
            taus[i] = T(0);
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::calcDetQ() const noexcept -> Tv {
        if constexpr (isComplex) {
            T x = 1;
            for (auto tau : taus)
                if (!tau.isZero()) [[likely]]
                    x *= -tau / tau.conjugate();
            return x;
        }
        else
            return std::ranges::count_if(taus, [](T x) { return !x.isZero(); }) % 2 == 0 ? 1 : -1;
    }
    /**
     * Decompose matrix like A = QDT(no pivoting), or A = QDTP(poviting), where D is diagonal
     */
    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::toQDT(VectorND<Tr>& diagD) noexcept {
        const size_t length = taus.getLength();
        assert(diagD.getLength() == length);
        for (size_t i = 0; i < length; ++i) {
            if (working(i, i).isSubNormal()) {
                diagD[i] = 1;
                continue;
            }
            if constexpr (isComplex)
                assert(working(i, i).imag().isZero() && "[Error]: Householder QR should have real diagonals");
            diagD[i] = working(i, i).real();
            working.row(i).tail(i) *= reciprocal(diagD[i]);
        }
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::toQDT() -> VectorND<Tr> {
        VectorND<Tr> vecD(taus.getLength());
        toQDT(vecD);
        return vecD;
    }

    template<Scalar T, bool Pivot>
    T QRDecomp<T, Pivot>::det() const noexcept {
        T result = calcDetQ() * getMatrixR().det();
        if constexpr (Pivot)
            result *= perm.det();
        return result;
    }

    template<Scalar T, bool Pivot>
    void QRDecomp<T, Pivot>::resize(size_t row, size_t col) {
        working.resize(row, col);
        auto l = std::min(row, col);
        taus.resize(l);
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
    auto QRDecomp<T, Pivot>::getMatrixQ() const -> MatrixND<T> {
        if constexpr (HasMKL()) {
            if constexpr (isComplex) // Our complex householder is slightly different from MKL, getMatrixQ_base() cannot apply to compute_mkl.
                return getMatrixQ_mkl();

            if (working.getSize() > Threhold)
                return getMatrixQ_mkl();
            return getMatrixQ_base();
        }
        else
            return getMatrixQ_base();
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::getMatrixQ_base() const -> MatrixND<T> {
        auto result = MatrixND<T>::unitMatrix(getRow());
        for (size_t i = 0; i < taus.getLength() - !working.isOverdetermined(); ++i) {
            auto block = result.rightCols(i);
            const auto col = working.col(i);
            applyHouseholder(block, taus[i], col.tail(i + 1));
        }
        return result;
    }

    template<Scalar T, bool Pivot>
    auto QRDecomp<T, Pivot>::getMatrixR() const noexcept {
        return working.triu();
    }
}

#ifdef PHYSICA_MKL
    #include "QRDecomp_MKL.h"
#endif
