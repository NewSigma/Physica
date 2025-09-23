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

#include "Physica/Core/Math/Statistics/NumCharacter.h"
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
        using MatrixND = Params::MatrixND;

        using Tr = T::RealType;
        using Trv = Tr::ValueType;
    private:
        const Params& params;
        DenseMatrix<Trv> aux;
        CyclicChainQDT<T> chainU;
        CyclicChainQDT<T> chainD;
        QRDecomp<T> kinetic;

        QRDecomp<T> qr;
        DiagMatrix<Tr> diagB;
        DiagMatrix<Tr> diagS;
        MatrixND buffer;

        Tr lnAbsDet;
        Tr signU = 1;
        Tr signD = 1;
        Array<MatrixND> greenUs;
        Array<MatrixND> greenDs;
    public:
        DQMC() = delete;
        DQMC(const Params& params_, int period);
        DQMC(const This&) = default;
        DQMC(This&&) noexcept = default;
        ~DQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        void step_random();

        template<RNG R>
        void step_mh();
        template<RNG R>
        void step_mh_for(int numStep);

        template<RNG R>
        [[nodiscard]] Trv step_spin();
        template<RNG R>
        void step_spin_for(int numStep);

        void invalidates();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params.getNumSplit(); }
        [[nodiscard]] auto& getAuxField() noexcept { return aux; }

        [[nodiscard]] Tr getLnAbsDet() const noexcept { return lnAbsDet; }
        [[nodiscard]] Tr getSignU() const noexcept { return signU; }
        [[nodiscard]] Tr getSignD() const noexcept { return signD; }
        [[nodiscard]] Tr getSign() const noexcept { return signU * signD; }
        [[nodiscard]] int getNumEqualGreen() const noexcept { return greenUs.getLength(); }
        [[nodiscard]] int getPeriod() const noexcept { return getNumSplit() / getNumEqualGreen(); }
        [[nodiscard]] const auto& getGreenUs() const noexcept { return greenUs; }
        [[nodiscard]] const auto& getGreenDs() const noexcept { return greenDs; }
    private:
        /* Operations */
        void resize(int numSite, int numSplit, int period);

        Vector2D<Tr> calcDelta(int site, int split) const noexcept;
        Vector2D<Tr> calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept;
        void single_flip(int site, int split) noexcept;
        void flipGreens(int site, int split, Vector2D<Tr> deltaRatios);
        void stepMHImpl(int site, int split, Tr prop);
        Trv stepSpinImpl(int site, int split, Tr prop);

        void splitDiag(const QDTDecomp<T>& qdt) noexcept;
        [[nodiscard]] std::array<Tr, 2> calcDet(const QDTDecomp<T>& qdt);
        [[nodiscard]] void calcGreens();
        /* Static members */
        template<RNG R>
        [[nodiscard]] static Array<int> makeRandomSites(int numSite);
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_, int period) : params(params_) {
        resize(params.getNumSite(), params.getNumSplit(), period);

        kinetic.compute(params.getExpB());
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
    void DQMC<T>::step_mh() {
        const int numSite = getNumSite();
        const auto props = VectorND<Tr>::template random_uniform<R>(numSite);
        const auto sites = makeRandomSites<R>(numSite);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(R::getInstance());
        for (int i = 0; i < numSite; ++i)
            stepMHImpl(sites[i], split, props[i]);

        calcGreens();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_mh_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        const int numSite = getNumSite();
        auto props = VectorND<Tr>(numSite);
        auto sites = makeRandomSites<R>(numSite);
        for (int step = 0; step < numStep; ++step) {
            props.template random_uniform<R>();
            int split = step % getNumSplit();
            for (int i = 0; i < numSite; ++i)
                stepMHImpl(sites[i], split, props[i]);

            calcGreens();
            std::ranges::shuffle(sites, R::getInstance());
        }
    }

    template<Scalar T>
    template<RNG R>
    auto DQMC<T>::step_spin() -> Trv {
        const int numSite = getNumSite();
        const auto props = VectorND<Tr>::template random_uniform<R>(numSite);
        const auto sites = makeRandomSites<R>(numSite);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(R::getInstance());
        Tr factor = 0;
        for (int i = 0; i < numSite; ++i)
            factor += stepSpinImpl(sites[i], split, props[i]);

        calcGreens();
        return factor;
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        const int numSite = getNumSite();
        auto props = VectorND<Tr>(numSite);
        auto sites = makeRandomSites<R>(numSite);
        for (int step = 0; step < numStep; ++step) {
            props.template random_uniform<R>();
            int split = step % getNumSplit();
            for (int i = 0; i < numSite; ++i)
                stepSpinImpl(sites[i], split, props[i]);

            calcGreens();
            std::ranges::shuffle(sites, R::getInstance());
        }
    }

    template<Scalar T>
    void DQMC<T>::invalidates() {
        const int numSplit = getNumSplit();
        DiagMatrix<Tr> expU(getNumSite());
        for (int split = 0; split < numSplit; ++split) {
            expU.diag() = exp(params.getAlpha() * aux.col(split));
            chainU[split].setMatrixR(kinetic.getMatrixR() * expU);
            expU.diag() = exp(-params.getAlpha() * aux.col(split));
            chainD[split].setMatrixR(kinetic.getMatrixR() * expU);
        }
        chainU.invalidates();
        chainD.invalidates();

        calcGreens();
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

        lnAbsDet.swap(obj.lnAbsDet);
        signU.swap(obj.signU);
        signD.swap(obj.signD);
        greenUs.swap(obj.greenUs);
        greenDs.swap(obj.greenDs);
    }

    template<Scalar T>
    void DQMC<T>::resize(int numSite, int numSplit, int period) {
        assert(numSite > 0 && "[Error]: Invalid NumSite");
        assert(numSplit > 0);
        assert((period > 0 && getNumSplit() % period == 0) && "[Error]: Invalid period");
        aux.resize(numSite, numSplit);
        chainU.resize(numSplit);
        chainD.resize(numSplit);
        kinetic.resize(numSite, numSite);

        qr.resize(numSite, numSite);
        diagB.resize(numSite);
        diagS.resize(numSite);
        buffer.resize(numSite, numSite);

        int numEqualGreen = numSplit / period;
        greenUs.resize(numEqualGreen, numSite, numSite);
        greenDs.resize(numEqualGreen, numSite, numSite);
    }

    template<Scalar T>
    auto DQMC<T>::calcDelta(int site, int split) const noexcept -> Vector2D<Tr> {
        const Tr x = Tr(2) * params.getAlpha() * aux(site, split);
        return exp(Vector2D<Tr>{-x, x}) - Tr(1); // The Eq. above Eq.(7.33) of [2]
    }

    template<Scalar T>
    auto DQMC<T>::calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept -> Vector2D<Tr> {
        assert(site < getNumSite() && split < getNumSplit());
        assert(getNumEqualGreen() == getNumSplit() && "[Error]: Cannot use rank-1 update if equal-time greens are not complete");
        const int split1 = (split + 1) % getNumSplit();
        Vector2D<Tr> result = deltas;
        result[0] *= Tr(1) - greenUs[split1](site, site).real();
        result[1] *= Tr(1) - greenDs[split1](site, site).real();
        result += Tr(1);
        return result; // Eq.(7.36) of [2]
    }

    template<Scalar T>
    void DQMC<T>::single_flip(int site, int split) noexcept {
        auto spins = aux.row(site);
        spins[split] = -spins[split];

        const Tr x = Tr(2) * params.getAlpha() * spins[split];
        const Vector2D<Tr> arr = exp(Vector2D<Tr>{x, -x});
        chainU.single_flip(site, split, arr[0], arr[1]);
        chainD.single_flip(site, split, arr[1], arr[0]);
    }

    template<Scalar T>
    void DQMC<T>::flipGreens(int site, int split, Vector2D<Tr> deltaRatios) {
        // Eq. (7.44) of [2]
        const int numSite = getNumSite();
        VectorND<T> vc(numSite);
        VectorND<T> vr(numSite);
        auto flipGreen = [site, &vc, &vr](MatrixND& green, Tr deltaRatio) {
            vc = (green - UnitMatrix<T>(green)).col(site);
            vr = green.row(site);
            green += deltaRatio * (vc * vr.transpose());
        };
        const int split1 = (split + 1) % getNumSplit();
        flipGreen(greenUs[split1], deltaRatios[0]);
        flipGreen(greenDs[split1], deltaRatios[1]);
    }

    template<Scalar T>
    void DQMC<T>::stepMHImpl(int site, int split, Tr prop) {
        const auto deltas = calcDelta(site, split);
        const auto ratios = calcRatio(site, split, deltas);
        if (prop < abs(ratios.prod())) {
            single_flip(site, split);
            flipGreens(site, split, divide(deltas, ratios));
        }
    }

    template<Scalar T>
    auto DQMC<T>::stepSpinImpl(int site, int split, Tr prop) -> Trv {
        const auto spins = aux.row(site);
        int spinUp = 0;
        for (int i = 0; i < spins.getLength(); ++i)
            spinUp += spins[i].isPositive();

        const auto& lnSpinWeights = params.getLnSpinWeights();
        const Trv lnWeight0 = lnSpinWeights[spinUp];
        const Trv lnWeight1 = lnSpinWeights[spinUp - aux(site, split).isPositive()];

        const auto deltas = calcDelta(site, split);
        const auto ratios = calcRatio(site, split, deltas);
        bool accept = prop < abs(ratios.prod()) * exp(lnWeight0 - lnWeight1);
        if (accept) {
            single_flip(site, split);
            flipGreens(site, split, divide(deltas, ratios));
        }
        return accept ? lnWeight1 : lnWeight0;
    }

    template<Scalar T>
    void DQMC<T>::splitDiag(const QDTDecomp<T>& qdt) noexcept {
        const auto& diagD = qdt.getMatrixD().diag();
        const Tr expBetaMu = exp(params.calcBetaMu());
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
    auto DQMC<T>::calcDet(const QDTDecomp<T>& qdt) -> std::array<Tr, 2> {
        splitDiag(qdt);

        buffer = qdt.getMatrixQ() * diagB.inverse();
        qr.compute(buffer.transpose() + diagS * qdt.getMatrixT());
        qr.getWorking().diag() += Tr(std::numeric_limits<T>::min()); // Handle potential underflow

        Tr lnAD = diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        Tr sign = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag().reals()).prod();
        assert(abs(sign) == Trv(1) && "[Error]: Bad sign");
        return {lnAD, sign};
    }

    template<Scalar T>
    void DQMC<T>::calcGreens() {
        const int numSite = getNumSite();
        auto temp = MatrixND(numSite, numSite);
        auto calcGreen = [this, &temp](const QDTDecomp<T>& qdt, MatrixND& green) -> std::array<Tr, 2> {
            const auto det = calcDet(qdt);
            temp = buffer * qr.getMatrixQ();
            green = qr.getMatrixR().inverse() * temp.transpose();
            return det;
        };

        const int numSplit = getNumSplit();
        const int period = getPeriod();
        Tr lnAbsDetU = 0;
        for (size_t i = 0; i < greenUs.getLength(); ++i) {
            const int shift = period * i;
            const int to = (numSplit + shift - 1) % numSplit;
            auto [lnAD, sign] = calcGreen(chainU.multiply(shift, to), greenUs[i]);

            assert((shift == 0 || signU == sign) && "[Error]: Unexpected sign mismatch");
            toNextMean(lnAbsDetU, i, lnAD);
            signU = sign;
        }

        Tr lnAbsDetD = 0;
        for (size_t i = 0; i < greenDs.getLength(); ++i) {
            const int shift = period * i;
            const int to = (numSplit + shift - 1) % numSplit;
            auto [lnAD, sign] = calcGreen(chainD.multiply(shift, to), greenDs[i]);

            assert((shift == 0 || signD == sign) && "[Error]: Unexpected sign mismatch");
            toNextMean(lnAbsDetD, i, lnAD);
            signD = sign;
        }
        lnAbsDet = lnAbsDetU + lnAbsDetD;
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
}
