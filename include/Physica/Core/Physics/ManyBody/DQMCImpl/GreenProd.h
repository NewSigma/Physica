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

#include "CyclicChainQDT.h"
#include "ImagKinetic.h"

namespace Physica {
    template<Scalar T>
    class GreenProd {
        using This = GreenProd<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
    private:
        MatrixND<T> expT;
        CyclicChainQDT<T> chainU;
        CyclicChainQDT<T> chainD;

        QRDecomp<T> qr;
        DiagMatrix<Tr> diagB;
        DiagMatrix<Tr> diagS;
        MatrixND<T>  buffer;

        Tr lnAbsDet;
        Tv sgnDet = 1;
    public:
        GreenProd() = delete;
        GreenProd(const HubbardParams<T>& params);
        GreenProd(const This&) = default;
        GreenProd(This&&) noexcept = default;
        ~GreenProd() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void invalidate(const MatrixND<Trv>& aux, Tr alpha, int split);
        void invalidates(const MatrixND<Trv>& aux, Tr alpha);

        void single_flip(int site, int split, Vector2D<Tr> factors, Vector2D<Tv> ratios) noexcept;
        void calcGreens(GreenPair& greens, int split, Tr betaMu);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return qr.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return chainU.getNumSplit(); }
    private:
        void splitDiag(const QDTDecomp<T>& qdt, Tr betaMu) noexcept;
        [[nodiscard]] std::pair<Tr, Tv> calcGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu);
    };

    template<Scalar T>
    GreenProd<T>::GreenProd(const HubbardParams<T>& params)
            : expT(params.getExpT())
            , chainU(params.getNumSplit())
            , chainD(params.getNumSplit())
            , qr(params.getNumSite(), params.getNumSite())
            , diagB(params.getNumSite())
            , diagS(params.getNumSite())
            , buffer(params.getNumSite()) {}

    template<Scalar T>
    void GreenProd<T>::invalidate(const MatrixND<Trv>& aux, Tr alpha, int split) {
        DiagMatrix<Tr> expU(getNumSite());
        expU.diag() = exp(alpha * aux.col(split));
        chainU[split] = expT * expU;
        expU.diag() = exp(-alpha * aux.col(split));
        chainD[split] = expT * expU;

        chainU.invalidate(split);
        chainD.invalidate(split);
    }

    template<Scalar T>
    void GreenProd<T>::invalidates(const MatrixND<Trv>& aux, Tr alpha) {
        const int numSplit = getNumSplit();
        DiagMatrix<Tr> expU(getNumSite());
        for (int split = 0; split < numSplit; ++split) {
            expU.diag() = exp(alpha * aux.col(split));
            chainU[split] = expT * expU;
            expU.diag() = exp(-alpha * aux.col(split));
            chainD[split] = expT * expU;
        }
        chainU.invalidates();
        chainD.invalidates();
    }

    template<Scalar T>
    void GreenProd<T>::single_flip(int site, int split, Vector2D<Tr> factors, Vector2D<Tv> ratios) noexcept {
        chainU.single_flip(site, split, factors[0], factors[1]);
        chainD.single_flip(site, split, factors[1], factors[0]);

        Tv prod = ratios.prod();
        lnAbsDet += ln(abs(prod));
        sgnDet *= unit(prod);
    }

    template<Scalar T>
    void GreenProd<T>::calcGreens(GreenPair& greens, int split, Tr betaMu) {
        const int numSplit = getNumSplit();
        const int from = (split + 1) % numSplit;
        const int to = (numSplit + split) % numSplit;
        Tr lnAbsDetU, lnAbsDetD;
        T signU, signD;
        {
            auto [lnAD, sign] = calcGreen(chainU.multiply(from, to), greens[0], betaMu);
            lnAbsDetU = lnAD;
            signU = sign;
        }
        {
            auto [lnAD, sign] = calcGreen(chainD.multiply(from, to), greens[1], betaMu);
            lnAbsDetD = lnAD;
            signD = sign;
        }
        lnAbsDet = lnAbsDetU + lnAbsDetD;
        sgnDet = signU * signD;
    }

    template<Scalar T>
    void GreenProd<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        chainU.swap(obj.chainU);
        chainD.swap(obj.chainD);

        qr.swap(obj.qr);
        diagB.swap(obj.diagB);
        diagS.swap(obj.diagS);
        buffer.swap(obj.buffer);
    }

    template<Scalar T>
    void GreenProd<T>::splitDiag(const QDTDecomp<T>& qdt, Tr betaMu) noexcept {
        const auto& diagD = qdt.getMatrixD().diag();
        const Tr expBetaMu = exp(betaMu);
        for (int i = 0; i < getNumSite(); ++i) {
            const Tr originD = diagD[i];
            const Tr expBetaMuD = expBetaMu * originD;
            const Tr absBetaMuD = abs(expBetaMuD);
            const bool sep = absBetaMuD > Tr(1);
            diagB.diag()[i] = sep ? absBetaMuD : Tr(1);
            diagS.diag()[i] = sep ? unit(originD) : expBetaMuD; // Use originD to avoid underflow
        }
    }

    template<Scalar T>
    auto GreenProd<T>::calcGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu) -> std::pair<Tr, Tv> {
        splitDiag(qdt, betaMu);

        buffer = qdt.getMatrixQ() * diagB.inv();
        qr.compute(buffer.hermite() + diagS * qdt.getMatrixT());
        qr.getWorking().diag() += Tr(std::numeric_limits<T>::min()); // Handle potential underflow

        Tr lnAD = diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        Tv sign = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag().reals()).prod();
        assert(T::isComplex || abs(sign) == Trv(1) && "[Error]: Bad sign");

        MatrixND<T> temp = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inv() * temp.hermite();
        return {lnAD, sign};
    }
}
