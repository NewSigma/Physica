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

#include "DQMCImpl/CyclicChainQDT.h"

namespace Physica {
    /**
     * Reference:
     * [1] Phys. Rev. B 40, 506; https://doi.org/10.1103/PhysRevB.40.506
     * [2] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016
     */
    template<Scalar T>
    class DQMC {
        using This = DQMC<T>;
        using Params = HubbardParams<T>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        constexpr static bool isComplex = T::isComplex;
    public:
        using MatrixND = Params::MatrixND;
        using GreenArray = Array2D<MatrixND, MatrixOption::Col, 2>;
    private:
        const Params* params;
        DenseMatrix<Trv> aux;
        CyclicChainQDT<T> chainU;
        CyclicChainQDT<T> chainD;
        QRDecomp<T> kinetic;

        QRDecomp<T> qr;
        DiagMatrix<Tr> diagB;
        DiagMatrix<Tr> diagS;
        MatrixND buffer;

        GreenArray greens;
        Tr lnAbsDet;
        Tv sgnDet = 1;
    public:
        DQMC() = delete;
        DQMC(const Params& params_);
        DQMC(const This&) = default;
        DQMC(This&&) noexcept = default;
        ~DQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        void step_random();

        template<RNG R>
        [[nodiscard]] const GreenArray& step_spin();
        template<RNG R>
        void step_spin_for(int numStep);

        template<RNG R>
        [[nodiscard]] std::pair<GreenArray, Vector2D<T>> step_pair();
        template<RNG R>
        void step_pair_for(int numStep);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return aux.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return aux.getCol(); }
        [[nodiscard]] const auto& getParams() const noexcept { return *params; }
        [[nodiscard]] const auto& getAuxField() const noexcept { return aux; }

        [[nodiscard]] Tr getLnAbsDet() const noexcept { return lnAbsDet; }
        [[nodiscard]] Tv getSign() const noexcept { return sgnDet; }
    private:
        /* Operations */
        void resize(int numSite, int numSplit);
        void invalidate(int split);
        void invalidates();

        Vector2D<Tr> calcDelta(int site, int split) const noexcept;
        Vector2D<Tv> calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept;
        void flipGreens(int site, int split, Vector2D<Tv> deltaRatios);
        void single_flip(int site, int split) noexcept;
        void single_flip(int site, int split, Vector2D<Tv> deltas, Vector2D<Tv> ratios) noexcept;
        void stepSpinImpl(int site, int split, Tr prob);
        void stepPairImpl(int pair, int site, int split, Tr prob);

        void splitDiag(const QDTDecomp<T>& qdt) noexcept;
        [[nodiscard]] std::pair<Tr, Tv> calcDet(const QDTDecomp<T>& qdt);
        void calcGreens();
        /* Static members */
        template<RNG R>
        [[nodiscard]] static Array<int> makeRandomSites(int numSite);
        static void checkSign(Tv sign1, Tv sign2) noexcept;
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_)
            : params(&params_)
            , chainU(params_.getNumSplit())
            , chainD(params_.getNumSplit()) {
        resize(params->getNumSite(), params->getNumSplit());

        kinetic.compute(params->getExpB());
        for (int split = 0; split < getNumSplit(); ++split) {
            chainU[split] = QDTDecomp<T>(kinetic);
            chainD[split] = QDTDecomp<T>(kinetic);
        }
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_random() {
        aux.template random_uniform<R>();
        aux = unit_elem(aux - Tr(0.5));
        invalidates();
    }

    template<Scalar T>
    template<RNG R>
    auto DQMC<T>::step_spin() -> const GreenArray& {
        const int numSite = getNumSite();
        const auto probs = VectorND<Tr>::template random_uniform<R>(numSite);
        const auto sites = makeRandomSites<R>(numSite);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(R::getInstance());
        for (int i = 0; i < numSite; ++i)
            stepSpinImpl(sites[i], split, probs[i]);

        calcGreens();
        return greens;
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        const int numSite = getNumSite();
        auto probs = VectorND<Tr>(numSite);
        auto sites = makeRandomSites<R>(numSite);
        auto dist = std::uniform_int_distribution<>(0, getNumSplit() - 1);
        for (int _ = 0; _ < numStep; ++_) {
            probs.template random_uniform<R>();
            int split = dist(R::getInstance());
            for (int i = 0; i < numSite; ++i)
                stepSpinImpl(sites[i], split, probs[i]);
            std::ranges::shuffle(sites, R::getInstance());
            calcGreens();
        }
    }

    template<Scalar T>
    template<RNG R>
    auto DQMC<T>::step_pair() -> std::pair<GreenArray, Vector2D<T>> {
        const int numSite = getNumSite();
        const auto probs = VectorND<Tr>::template random_uniform<R>(numSite);
        const auto sites = makeRandomSites<R>(numSite);
        const int pair = std::uniform_int_distribution<>(0, getNumSite() - 1)(R::getInstance());
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(R::getInstance());
        for (int i = 0; i < numSite; ++i) {
            int site = sites[i];
            if (pair == site) [[unlikely]]
                continue;
            stepPairImpl(pair, site, split, probs[i]);
        }

        const int split1 = (split + 1) % getNumSplit();
        GreenArray result(2, 2);
        for (int i : {0, 1})
            result(i, 0) = greens(i, split1);

        auto weights = [this, pair, split]() -> Vector2D<T> {
            const auto deltas = calcDelta(pair, split);
            const auto ratios = calcRatio(pair, split, deltas);
            auto w0 = ratios.prod();
            auto w1 = Trv(1) + w0;
            auto absw = abs(w1);
            auto sign = getSign();
            lnAbsDet += ln(absw);
            sgnDet = sign * unit(w1);

            single_flip(pair, split);
            return {sign / absw, sign * w0 / absw};
        }();

        for (int i : {0, 1})
            result(i, 1) = greens(i, split1);
        calcGreens();
        return std::make_pair(std::move(result), std::move(weights));
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_pair_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        const int numSite = getNumSite();
        auto probs = VectorND<Tr>(numSite);
        auto sites = makeRandomSites<R>(numSite);
        auto dist0 = std::uniform_int_distribution<>(0, getNumSite() - 1);
        auto dist1 = std::uniform_int_distribution<>(0, getNumSplit() - 1);
        for (int _ = 0; _ < numStep; ++_) {
            probs.template random_uniform<R>();
            int pair = dist0(R::getInstance());
            int split = dist1(R::getInstance());
            for (int i = 0; i < numSite; ++i) {
                int site = sites[i];
                if (pair == site) [[unlikely]]
                    continue;
                stepPairImpl(pair, site, split, probs[i]);
            }
            std::ranges::shuffle(sites, R::getInstance());
            calcGreens();
        }
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        aux.swap(obj.aux);
        chainU.swap(obj.chainU);
        chainD.swap(obj.chainD);
        kinetic.swap(obj.kinetic);

        qr.swap(obj.qr);
        diagB.swap(obj.diagB);
        diagS.swap(obj.diagS);
        buffer.swap(obj.buffer);

        greens.swap(obj.greens);
        lnAbsDet.swap(obj.lnAbsDet);
        sgnDet.swap(obj.sgnDet);
    }

    template<Scalar T>
    void DQMC<T>::resize(int numSite, int numSplit) {
        assert(numSite > 0 && "[Error]: Invalid NumSite");
        assert(numSplit > 0);
        aux.resize(numSite, numSplit);
        kinetic.resize(numSite, numSite);

        qr.resize(numSite, numSite);
        diagB.resize(numSite);
        diagS.resize(numSite);
        buffer.resize(numSite, numSite);

        greens.resize(2, numSplit, numSite, numSite);
    }

    template<Scalar T>
    void DQMC<T>::invalidate(int split) {
        DiagMatrix<Tr> expU(getNumSite());
        expU.diag() = exp(params->getAlpha() * aux.col(split));
        chainU[split].setMatrixR(kinetic.getMatrixR() * expU);
        expU.diag() = exp(-params->getAlpha() * aux.col(split));
        chainD[split].setMatrixR(kinetic.getMatrixR() * expU);

        chainU.invalidate(split);
        chainD.invalidate(split);

        calcGreens();
    }

    template<Scalar T>
    void DQMC<T>::invalidates() {
        const int numSplit = getNumSplit();
        DiagMatrix<Tr> expU(getNumSite());
        for (int split = 0; split < numSplit; ++split) {
            expU.diag() = exp(params->getAlpha() * aux.col(split));
            chainU[split].setMatrixR(kinetic.getMatrixR() * expU);
            expU.diag() = exp(-params->getAlpha() * aux.col(split));
            chainD[split].setMatrixR(kinetic.getMatrixR() * expU);
        }
        chainU.invalidates();
        chainD.invalidates();

        calcGreens();
    }

    template<Scalar T>
    auto DQMC<T>::calcDelta(int site, int split) const noexcept -> Vector2D<Tr> {
        const Tr x = Trv(2) * params->getAlpha() * aux(site, split);
        return exp(Vector2D<Tr>{-x, x}) - Trv(1); // The Eq. above Eq.(7.33) of [2]
    }

    template<Scalar T>
    auto DQMC<T>::calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept -> Vector2D<Tv> {
        assert(site < getNumSite() && split < getNumSplit());
        const int split1 = (split + 1) % getNumSplit();
        Vector2D<Tv> result = deltas;
        result[0] *= Trv(1) - greens(0, split1)(site, site);
        result[1] *= Trv(1) - greens(1, split1)(site, site);
        result += Trv(1);
        return result; // Eq.(7.36) of [2]
    }

    template<Scalar T>
    void DQMC<T>::flipGreens(int site, int split, Vector2D<Tv> deltaRatios) {
        // Eq. (7.44) of [2]
        const int numSite = getNumSite();
        VectorND<T> vc(numSite);
        VectorND<T> vr(numSite);
        auto flipGreen = [site, &vc, &vr](MatrixND& green, Tv deltaRatio) {
            vc = (green - UnitMatrix<T>(green)).col(site);
            vr = green.row(site);
            green += deltaRatio * (vc * vr.transpose());
        };
        const int split1 = (split + 1) % getNumSplit();
        flipGreen(greens(0, split1), deltaRatios[0]);
        flipGreen(greens(1, split1), deltaRatios[1]);
    }

    template<Scalar T>
    void DQMC<T>::single_flip(int site, int split) noexcept {
        const auto deltas = calcDelta(site, split);
        const auto ratios = calcRatio(site, split, deltas);
        single_flip(site, split, deltas, ratios);
    }

    template<Scalar T>
    void DQMC<T>::single_flip(int site, int split, Vector2D<Tv> deltas, Vector2D<Tv> ratios) noexcept {
        assert(deltas == calcDelta(site, split));
        assert(ratios == calcRatio(site, split, deltas));
        auto spins = aux.row(site);
        spins[split] = -spins[split];

        const Tr x = Tr(2) * params->getAlpha() * spins[split];
        const Vector2D<Tr> arr = exp(Vector2D<Tr>{x, -x});
        chainU.single_flip(site, split, arr[0], arr[1]);
        chainD.single_flip(site, split, arr[1], arr[0]);

        flipGreens(site, split, divide(deltas, ratios));

        Tv prod = ratios.prod();
        lnAbsDet += ln(abs(prod));
        sgnDet *= unit(prod);
    }

    template<Scalar T>
    void DQMC<T>::stepSpinImpl(int site, int split, Tr prob) {
        const auto deltas = calcDelta(site, split);
        const auto ratios = calcRatio(site, split, deltas);
        const bool accept = prob < abs(ratios.prod());
        if (accept)
            single_flip(site, split, deltas, ratios);
    }

    template<Scalar T>
    void DQMC<T>::stepPairImpl(int pair, int site, int split, Tr prob) {
        assert(pair != site);
        assert(pair < getNumSite() && site < getNumSite());
        const auto prob1 = [this, pair, site, split]() {
            const auto siteD = calcDelta(site, split);
            const auto siteR = calcRatio(site, split, siteD);

            const auto pairD = calcDelta(pair, split);
            const auto pairR0 = calcRatio(pair, split, pairD).prod();

            single_flip(site, split, siteD, siteR);
            const auto pairR1 = calcRatio(pair, split, pairD).prod();
            return abs(siteR.prod() * ((Trv(1) + pairR1) / (Trv(1) + pairR0)));
        }();

        if (!(prob < prob1))
            single_flip(site, split);
    }

    template<Scalar T>
    void DQMC<T>::splitDiag(const QDTDecomp<T>& qdt) noexcept {
        const auto& diagD = qdt.getMatrixD().diag();
        const Tr expBetaMu = exp(params->calcBetaMu());
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
    auto DQMC<T>::calcDet(const QDTDecomp<T>& qdt) -> std::pair<Tr, Tv> {
        splitDiag(qdt);

        buffer = qdt.getMatrixQ() * diagB.inv();
        qr.compute(buffer.hermite() + diagS * qdt.getMatrixT());
        qr.getWorking().diag() += Tr(std::numeric_limits<T>::min()); // Handle potential underflow

        Tr lnAD = diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        Tv sign = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag().reals()).prod();
        assert(isComplex || abs(sign) == Trv(1) && "[Error]: Bad sign");
        return {lnAD, sign};
    }
    // TODO: Add exception
    template<Scalar T>
    void DQMC<T>::calcGreens() {
        const int numSite = getNumSite();
        auto temp = MatrixND(numSite, numSite);
        auto calcGreen = [this, &temp](const QDTDecomp<T>& qdt, MatrixND& green) -> std::pair<Tr, Tv> {
            const auto det = calcDet(qdt);
            temp = buffer * qr.getMatrixQ();
            green = qr.getMatrixR().inv() * temp.hermite();
            return det;
        };

        const int numSplit = getNumSplit();
        Tr lnAbsDetU = 0;
        Tv signU;
        for (size_t i = 0; i < numSplit; ++i) {
            const int to = (numSplit + i - 1) % numSplit;
            auto [lnAD, sign] = calcGreen(chainU.multiply(i, to), greens(0, i));
            if (i != 0)
                checkSign(signU, sign);

            lnAbsDetU.toNextMean(i, lnAD);
            signU = sign;
        }

        Tr lnAbsDetD = 0;
        Tv signD;
        for (size_t i = 0; i < numSplit; ++i) {
            const int to = (numSplit + i - 1) % numSplit;
            auto [lnAD, sign] = calcGreen(chainD.multiply(i, to), greens(1, i));
            if (i != 0)
                checkSign(signD, sign);

            lnAbsDetD.toNextMean(i, lnAD);
            signD = sign;
        }
        lnAbsDet = lnAbsDetU + lnAbsDetD;
        sgnDet = signU * signD;
    }

    template<Scalar T>
    template<RNG R>
    Array<int> DQMC<T>::makeRandomSites(int numSite) {
        auto sites = Array<int>(numSite);
        for (int i = 0; i < numSite; ++i)
            sites[i] = i;
        std::ranges::shuffle(sites, R::getInstance());
        return sites;
    }

    template<Scalar T>
    void DQMC<T>::checkSign([[maybe_unused]] Tv sign1, [[maybe_unused]] Tv sign2) noexcept {
        // TODO: check complex sign
        if constexpr (!isComplex)
            assert(sign1 == sign2 && "[Error]: Unexpected sign mismatch");
    }
}
