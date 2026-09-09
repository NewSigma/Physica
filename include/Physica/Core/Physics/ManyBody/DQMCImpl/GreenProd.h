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
#include "Physica/Core/Parallel/Parallel.h"

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
        Array<CyclicChainQDT<T>, 2> chains;
        Array<DenseQR<T>, 2> qrs;
        Array<DiagMatrix<Tr>, 2> diagBs;
        Array<DiagMatrix<Tr>, 2> diagSs;
        Array<MatrixND<T>, 2> buffers;
    public:
        GreenProd() = delete;
        GreenProd(const HubbardParams<T>& params);
        GreenProd(const This&) = default;
        GreenProd(This&&) noexcept = default;
        ~GreenProd() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void invalidate(int split);
        void invalidates(const MatrixND<Trv>& aux, Tr alpha);

        void single_flip(int site, int split, Vector2D<Tr> factors) noexcept;
        template<ExecutePolicy P = Sequential>
        auto calcGreens(GreenPair& greens, int split, Tr betaMu);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return qrs.front().getOrder(); }
        [[nodiscard]] int getNumSplit() const noexcept { return chains.front().getNumSplit(); }
    private:
        void splitDiag(const QDTDecomp<T>& qdt, Tr betaMu, int spin) noexcept;
        void calcGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu, int spin);
        [[nodiscard]] std::pair<Tr, Tv> calcDetGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu, int spin);
    };

    template<Scalar T>
    GreenProd<T>::GreenProd(const HubbardParams<T>& params)
            : expT(params.getExpT())
            , chains(2, params.getNumSplit())
            , qrs(2, params.getNumSite(), params.getNumSite())
            , diagBs(2, params.getNumSite())
            , diagSs(2, params.getNumSite())
            , buffers(2, params.getNumSite()) {}

    template<Scalar T>
    void GreenProd<T>::invalidate(int split) {
        for (auto& chain : chains)
            chain.invalidate(split);
    }

    template<Scalar T>
    void GreenProd<T>::invalidates(const MatrixND<Trv>& aux, Tr alpha) {
        for (auto& chain : chains)
            chain.invalidates();

        const int numSplit = getNumSplit();
        DiagMatrix<Tr> expU(getNumSite());
        for (int split = 0; split < numSplit; ++split) {
            expU.diag() = exp(alpha * aux.col(split));
            chains[0][split] = expT * expU;
            expU.diag() = exp(-alpha * aux.col(split));
            chains[1][split] = expT * expU;
        }
    }

    template<Scalar T>
    void GreenProd<T>::single_flip(int site, int split, Vector2D<Tr> factors) noexcept {
        chains[0].single_flip(site, split, factors[0], factors[1]);
        chains[1].single_flip(site, split, factors[1], factors[0]);
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto GreenProd<T>::calcGreens(GreenPair& greens, int split, Tr betaMu) {
        const int numSplit = getNumSplit();
        const int from = (split + 1) % numSplit;
        const int to = (numSplit + split) % numSplit;
        return parallel_for<P>([this, &greens, from, to, betaMu](int spin) {
            const auto& qdt = chains[spin].multiply(from, to);
            calcGreen(qdt, greens[spin], betaMu, spin);
        }, 2);
    }

    template<Scalar T>
    void GreenProd<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        chains.swap(obj.chains);
        qrs.swap(obj.qrs);
        diagBs.swap(obj.diagBs);
        diagSs.swap(obj.diagSs);
        buffers.swap(obj.buffers);
    }

    template<Scalar T>
    void GreenProd<T>::splitDiag(const QDTDecomp<T>& qdt, Tr betaMu, int spin) noexcept {
        const auto& diagD = qdt.getMatrixD().diag();
        const Tr expBetaMu = exp(betaMu);
        auto& diagB = diagBs[spin];
        auto& diagS = diagSs[spin];
        for (int i = 0; i < getNumSite(); ++i) {
            const Tr expBetaMuD = expBetaMu * diagD[i];
            const Tr absBetaMuD = abs(expBetaMuD);
            const Tr halfAbs = sqrt(absBetaMuD);
            bool isSubNormal = halfAbs.isSubNormal();
            diagB.diag()[i] = isSubNormal ? Tr(1) : halfAbs;
            diagS.diag()[i] = isSubNormal ? Tr(0) : expBetaMuD / halfAbs;
        }
    }

    template<Scalar T>
    void GreenProd<T>::calcGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu, int spin) {
        splitDiag(qdt, betaMu, spin);

        auto& diagB = diagBs[spin];
        auto& diagS = diagSs[spin];
        auto& buffer = buffers[spin];
        auto& qr = qrs[spin];
        buffer = qdt.getMatrixQ() * diagB.inv();
        qr.compute(buffer.hermite() + diagS * qdt.getMatrixT());
        qr.getWorking().diag() += Tr(std::numeric_limits<T>::min()); // Handle potential underflow

        MatrixND<T> temp = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inv() * temp.hermite();
    }

    template<Scalar T>
    auto GreenProd<T>::calcDetGreen(const QDTDecomp<T>& qdt, MatrixND<T>& green, Tr betaMu, int spin) -> std::pair<Tr, Tv> {
        calcGreen(qdt, green, betaMu, spin);

        auto& diagB = diagBs[spin];
        auto& qr = qrs[spin];
        Tr lnAD = diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        Tv sgnD = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag().reals()).prod();
        assert(T::isComplex() || abs(sgnD) == Trv(1) && "[Error]: Bad sign");
        return {lnAD, sgnD};
    }
}
