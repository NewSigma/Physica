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
#include "DenseQR.h"

namespace Physica {
    template<Scalar Tv, bool Pivot>
    class DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot> {
        static_assert(!Tv::isDiffable(), "[Error]: Invalid tparam");
        static_assert(!Pivot, "[Error]: Pivot is not supported for forward diff");
        using T = Diff<Tv, DiffMode::Forward, 1>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;

        using This = DenseQR<T, Pivot>;
        constexpr static bool isComplex = T::isComplex();
    private:
        MatrixND<T> matrixQ;
        MatrixND<T> working;
        VectorND<Tv> taus;
        Tv detQ{};
    public:
        DenseQR() = default;
        explicit DenseQR(size_t order);
        DenseQR(size_t row, size_t col);
        DenseQR(const Matrix auto& source);
        DenseQR(const This&) = default;
        DenseQR(This&&) noexcept = default;
        ~DenseQR() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        template<Matrix M>
        void compute(const M& source);

        [[nodiscard]] Tv calcDetQ() const noexcept;
        void toQDT(VectorND<Tr>& diagD) noexcept;
        [[nodiscard]] VectorND<Tr> toQDT();

        [[nodiscard]] T det() const noexcept;

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getWorking(this auto&& self) noexcept;
        [[nodiscard]] const auto& getTaus() const noexcept { return taus; }
        [[nodiscard]] size_t getRow() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return working.getCol(); }
        [[nodiscard]] size_t getOrder() const noexcept { return working.getOrder(); }
        [[nodiscard]] MatrixND<T> getMatrixQ() const;
        [[nodiscard]] auto getMatrixR() const noexcept;
    };

    template<Scalar Tv, bool Pivot>
    DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::DenseQR(size_t order) : DenseQR(order, order) {}

    template<Scalar Tv, bool Pivot>
    DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::DenseQR(size_t row, size_t col) {
        resize(row, col);
    }

    template<Scalar Tv, bool Pivot>
    DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::DenseQR(const Matrix auto& source)
            : DenseQR(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar Tv, bool Pivot>
    template<Matrix M>
    void DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::compute(const M& source) {
        static_assert(std::is_same_v<typename M::ScalarType, T>, "[Error]: Inconsistent ScalarType");
        assert(getOrder() == source.getOrder());
        assert(source.isSquare() && "[Error]: Forward diff of QR requires a square matrix");
        const size_t order = getRow();
        {
            DenseQR<Tv, false> solver(order);
            solver.compute(source.values());
            resize(order, order);

            matrixQ.values() = std::move(solver.getMatrixQ());
            working.values() = solver.getMatrixR();
            taus = std::move(solver.getTaus());
            detQ = solver.calcDetQ();
        }

        const auto& q = matrixQ.values();
        const auto& r = working.values();
        const MatrixND<Tv> temp = q.hermite() * source.grads();
        const MatrixND<Tv> temp1 = temp * r.triu().inv();
        // TODO: Add triu(int offset); Add AntiHermiteMatrix
        const MatrixND<Tv> lower = temp1 - temp1.triu();
        const MatrixND<Tv> omega = lower - lower.hermite();
        matrixQ.grads() = q * omega;
        working.grads() = (temp - omega * r).triu();
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::calcDetQ() const noexcept -> Tv {
        return detQ;
    }

    template<Scalar Tv, bool Pivot>
    void DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::toQDT(VectorND<Tr>& diagD) noexcept {
        const size_t length = taus.getLength();
        assert(diagD.getLength() == length);
        for (size_t i = 0; i < length; ++i) {
            if (working[i, i].isSubNormal()) {
                diagD[i] = 1;
                continue;
            }
            if constexpr (isComplex) {
                [[maybe_unused]] bool isReal = abs(working[i, i].imag()) < abs(working[i, i].real()) * sqrt(Trv(std::numeric_limits<T>::epsilon()));
                assert(isReal && "[Error]: Householder QR should have real diagonals");
            }
            diagD[i] = working[i, i].real();
            working.row(i).tail(i) *= reciprocal(diagD[i]);
        }
        working.grads().diag().zeros();
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::toQDT() -> VectorND<Tr> {
        VectorND<Tr> vecD(taus.getLength());
        toQDT(vecD);
        return vecD;
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::det() const noexcept -> T {
        return detQ * getMatrixR().det();
    }

    template<Scalar Tv, bool Pivot>
    void DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::resize(size_t row, size_t col) {
        assert(row == col && "[Error]: Forward diff of QR requires a square matrix");
        matrixQ.resize(row, row);
        working.resize(row, col);
        taus.resize(std::min(row, col));
    }

    template<Scalar Tv, bool Pivot>
    void DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        matrixQ.swap(obj.matrixQ);
        working.swap(obj.working);
        taus.swap(obj.taus);
        detQ.swap(obj.detQ);
    }

    template<Scalar Tv, bool Pivot>
    auto&& DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::getWorking(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).working;
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::getMatrixQ() const -> MatrixND<T> {
        return matrixQ;
    }

    template<Scalar Tv, bool Pivot>
    auto DenseQR<Diff<Tv, DiffMode::Forward, 1>, Pivot>::getMatrixR() const noexcept {
        return working.triu();
    }
}
