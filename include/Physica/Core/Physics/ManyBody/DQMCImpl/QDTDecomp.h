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

#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/QRDecomp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "HubbardParams.h"

namespace Physica {
    template<Scalar T>
    class QDTDecomp {
        using This = QDTDecomp<T>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        MatrixND<T> matrixQ;
        DiagMatrix<Tr> matrixD; // matrixD suffers from over/underflow
        QRDecomp<T> qr;
        Tv detQ;
    public:
        QDTDecomp() = default;
        QDTDecomp(size_t size);
        QDTDecomp(QRDecomp<T> qr_);
        QDTDecomp(const Matrix auto& source);
        QDTDecomp(const This&) = default;
        QDTDecomp(This&&) noexcept = default;
        ~QDTDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] This operator*(const This& other) const;
        /* Operations */
        void compute(const Matrix auto& source);
        [[nodiscard]] Tv calcDetQ() const noexcept;

        void single_flip(int site, Tr factor, Tr invfac) noexcept;

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        const This& dump() const;
        /* Getters */
        [[nodiscard]] const auto& getMatrixQ() const noexcept { return matrixQ; }
        [[nodiscard]] const auto& getMatrixD() const noexcept { return matrixD; }
        [[nodiscard]] const auto& getMatrixT() const noexcept { return qr.getWorking(); }
        [[nodiscard]] auto& getMatrixT() noexcept { return qr.getWorking(); }
        [[nodiscard]] const auto& getQR() const noexcept { return qr; }
        [[nodiscard]] size_t getSize() const noexcept { return matrixQ.getRow(); }
        /* Setters */
        void setMatrixR(const Matrix auto& matrixR);
    private:
        void toQDT();
    };

    template<Scalar T>
    QDTDecomp<T>::QDTDecomp(size_t size) : matrixQ(size, size), matrixD(size), qr(size, size) {}

    template<Scalar T>
    QDTDecomp<T>::QDTDecomp(QRDecomp<T> qr_) : QDTDecomp(qr_.getRow()) {
        assert(qr_.getWorking().isSquare());
        qr = std::move(qr_);
        detQ = qr.calcDetQ();
        toQDT();
    }

    template<Scalar T>
    QDTDecomp<T>::QDTDecomp(const Matrix auto& source) : QDTDecomp(source.getRow()) {
        compute(source);
    }
    /**
     * Reference:
     * [1] Linear Algebra and its Applications 435(3), 659-673 (2011); https://doi.org/10.1016/j.laa.2010.06.023
     */
    template<Scalar T>
    auto QDTDecomp<T>::operator*(const This& other) const -> This {
        MatrixND<T> buffer = getMatrixT() * other.getMatrixQ();
        buffer = getMatrixD() * buffer;
        buffer = buffer * other.getMatrixD();

        auto result = This(buffer);
        buffer = getMatrixQ() * result.getMatrixQ();
        buffer.swap(result.matrixQ);

        buffer = result.getMatrixT() * other.getMatrixT();
        buffer.swap(result.qr.getWorking());

        result.detQ = detQ * result.qr.calcDetQ();
        return result;
    }

    template<Scalar T>
    void QDTDecomp<T>::compute(const Matrix auto& source) {
        assert(source.isSquare());
        if (source.getRow() != getSize()) [[unlikely]]
            resize(source.getRow());

        qr.compute(source);
        detQ = qr.calcDetQ();
        toQDT();
    }

    template<Scalar T>
    auto QDTDecomp<T>::calcDetQ() const noexcept -> Tv {
        return detQ;
    }

    template<Scalar T>
    void QDTDecomp<T>::single_flip(int site, Tr factor, Tr invfac) noexcept {
        assert(scalarNear(factor * invfac, Tr(1), std::numeric_limits<T>::epsilon() * 10) && "[Error]: Invalid argument");
        if (site > 0) [[likely]]
            getMatrixT().col(site).head(site) *= factor;
        if (site + 1 < getSize()) [[likely]]
            getMatrixT().row(site).tail(site + 1) *= invfac;
        matrixD.diag()[site] *= factor;
    }

    template<Scalar T>
    void QDTDecomp<T>::resize(size_t size) {
        matrixQ.resize(size, size);
        matrixD.resize(size);
        qr.resize(size, size);
    }

    template<Scalar T>
    void QDTDecomp<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        matrixQ.swap(obj.matrixQ);
        matrixD.swap(obj.matrixD);
        qr.swap(obj.qr);
        detQ.swap(obj.detQ);
    }

    template<Scalar T>
    auto QDTDecomp<T>::dump() const -> const This& {
        MatrixND<T> x = getMatrixQ() * getMatrixD() * getMatrixT();
        std::cout << x;
        return *this;
    }

    template<Scalar T>
    void QDTDecomp<T>::setMatrixR(const Matrix auto& matrixR) {
        qr.getWorking() = matrixR; // Clear lower diagonal
        qr.toQDT(matrixD.diag());
    }

    template<Scalar T>
    void QDTDecomp<T>::toQDT() {
        matrixQ = qr.getMatrixQ();
        setMatrixR(qr.getMatrixR());
    }
}
