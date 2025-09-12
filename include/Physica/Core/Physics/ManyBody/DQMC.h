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

        Tr lnPartitionZ;
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
        void step_spin();
        template<RNG R>
        void step_spin_for(int numStep);

        void invalidates();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params.getNumSplit(); }
        [[nodiscard]] auto& getAuxField() noexcept { return aux; }

        [[nodiscard]] Tr getLnPartitionZ() const noexcept { return lnPartitionZ; }
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

        void single_flip(int site, int split) noexcept;
        void splitDiag(const QDTDecomp<T>& qdt) noexcept;
        [[nodiscard]] std::pair<Tr, Tr> lnPartition(const QDTDecomp<T>& qdt);
        [[nodiscard]] std::pair<Tr, Tr> calcGreen(const QDTDecomp<T>& qdt, MatrixND& green);
        [[nodiscard]] Tr calcGreens();

        Vector2D<Tr> calcRatio(int site, int split) const noexcept;
        Tr calcLnSpinWaveWeight(int site) const noexcept;
        Tr calcLnSpinWaveWeight(Tr sumSpin) const noexcept;
        /* Static members */
        template<RNG R>
        static bool accept(Tr deltaW) noexcept;
        static bool accept(Tr deltaW, Tr prop) noexcept;
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
        auto& rng = R::getInstance();
        auto props = VectorND<Tr>(getNumSite());
        auto sites = Array<int>(getNumSite());
        for (int i = 0; i < getNumSite(); ++i)
            sites[i] = i;

        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(rng);
        std::shuffle(sites.begin(), sites.end(), rng);
        props.template random_uniform<R>();
        for (int i = 0; i < getNumSite(); ++i) {
            int site = sites[i];
            if (props[i] < abs(calcRatio(site, split).prod()))
                single_flip(site, split);
        }
        lnPartitionZ = calcGreens();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_mh_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        auto props = VectorND<Tr>(getNumSite());
        auto sites = Array<int>(getNumSite());
        for (int i = 0; i < getNumSite(); ++i)
            sites[i] = i;

        for (int step = 0; step < numStep; ++step) {
            std::shuffle(sites.begin(), sites.end(), R::getInstance());
            props.template random_uniform<R>();
            int split = step % getNumSplit();
            for (int i = 0; i < getNumSite(); ++i) {
                int site = sites[i];
                if (props[i] < abs(calcRatio(site, split).prod()))
                    single_flip(site, split);
            }
            lnPartitionZ = calcGreens();
        }
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin() {
        auto& rng = R::getInstance();
        const int site = std::uniform_int_distribution<>(0, getNumSite() - 1)(rng);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(rng);
        auto spins = aux.row(site);
        const Tr sumSpin0 = spins.sum();
        const Tr sumSpin1 = sumSpin0 - Tr(2) * spins[split];
        const Tr deltaW = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
        const bool flip1 = accept<R>(deltaW);
        if (flip1)
            single_flip(site, split);

        const auto p = Tr::template random_uniform<R>();
        const bool flip2 = p < 0.5;
        if (flip2) {
            for (int i = 0; i < getNumSplit(); ++i)
                single_flip(site, i);
        }

        const bool noFlip = !flip1 && !flip2;
        if (noFlip)
            return;
    
        lnPartitionZ = calcGreens();
        for (int i = 0; i < getNumSite(); ++i)
            lnPartitionZ -= calcLnSpinWaveWeight(i);
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin_for(int numStep)  {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        auto props = VectorND<Tr>(getNumSplit());
        auto splits = Array<int>(getNumSplit());
        for (int i = 0; i < getNumSplit(); ++i)
            splits[i] = i;

        for (int site = 0; site < getNumSite(); ++site) {
            auto spins = aux.row(site);
            Tr sumSpin0 = spins.sum();
            Tr weight0 = calcLnSpinWaveWeight(sumSpin0);
            for (int step = 0; step < numStep; ++step) {
                const int i = step % getNumSplit();
                if (i == 0) {
                    props.template random_uniform<R>();
                    std::shuffle(splits.begin(), splits.end(), R::getInstance());
                }

                const int split = splits[i];
                const Tr sumSpin1 = sumSpin0 - Tr(2) * spins[split];
                const Tr weight1 = calcLnSpinWaveWeight(sumSpin1);
                const Tr deltaW = weight1 - weight0;
                if (accept(deltaW, props[i])) {
                    spins[split] = -spins[split];
                    sumSpin0 = sumSpin1;
                    weight0 = weight1;
                }
            }

            const auto p = Tr::template random_uniform<R>();
            if (p < 0.5)
                spins = -spins;
        }

        invalidates();
        for (int i = 0; i < getNumSite(); ++i)
            lnPartitionZ -= calcLnSpinWaveWeight(i);
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

        lnPartitionZ = calcGreens();
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

        lnPartitionZ.swap(obj.lnPartitionZ);
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
    void DQMC<T>::single_flip(int site, int split) noexcept {
        auto spins = aux.row(site);
        spins[split] = -spins[split];
        // TODO: Make it a data member
        std::array<Tr, 2> arr{exp(Trv(2) * params.getAlpha() * spins[split]), exp(Trv(-2) * params.getAlpha() * spins[split])};
        chainU.single_flip(site, split, arr[0], arr[1]);
        chainD.single_flip(site, split, arr[1], arr[0]);
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
    auto DQMC<T>::lnPartition(const QDTDecomp<T>& qdt) -> std::pair<Tr, Tr> {
        splitDiag(qdt);

        buffer = qdt.getMatrixQ() * diagB.inverse();
        qr.compute(buffer.transpose() + diagS * qdt.getMatrixT());
        qr.getWorking().diag() += Tr(std::numeric_limits<T>::min()); // Handle potential underflow

        Tr lnZ = diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        Tr sign = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag().reals()).prod();
        assert(abs(sign) == Trv(1) && "[Error]: Bad sign");
        return {lnZ, sign};
    }

    template<Scalar T>
    auto DQMC<T>::calcGreen(const QDTDecomp<T>& qdt, MatrixND& green) -> std::pair<Tr, Tr> {
        const auto pair = lnPartition(qdt);
        MatrixND temp = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inverse() * temp.transpose();
        return pair;
    }

    template<Scalar T>
    auto DQMC<T>::calcGreens() -> Tr {
        const int numSplit = getNumSplit();
        const int period = getPeriod();
        Tr lnZU = 0;
        for (size_t i = 0; i < greenUs.getLength(); ++i) {
            const int shift = period * i;
            const int to = (numSplit + shift - 1) % numSplit;
            auto [lnZ, sign] = calcGreen(chainU.multiply(shift, to), greenUs[i]);

            assert((shift == 0 || signU == sign) && "[Error]: Unexpected sign mismatch");
            toNextMean(lnZU, i, lnZ);
            signU = sign;
        }

        Tr lnZD = 0;
        for (size_t i = 0; i < greenDs.getLength(); ++i) {
            const int shift = period * i;
            const int to = (numSplit + shift - 1) % numSplit;
            auto [lnZ, sign] = calcGreen(chainD.multiply(shift, to), greenDs[i]);

            assert((shift == 0 || signD == sign) && "[Error]: Unexpected sign mismatch");
            toNextMean(lnPartitionZ, i, lnZ);
            signD = sign;
        }
        return lnZU + lnZD;
    }

    template<Scalar T>
    auto DQMC<T>::calcRatio(int site, int split) const noexcept -> Vector2D<Tr> {
        assert(site < getNumSite() && split < getNumSplit());
        assert(getNumEqualGreen() == getNumSplit() && "[Error]: Cannot use rank-1 update if equal-time greens are not complete");
        const Tr delta = Tr(2) * params.getAlpha() * aux(site, split);
        Vector2D<Tr> result{-delta, delta};
        result = exp(result) - Tr(1);
        result[0] *= Tr(1) - greenUs[split](site, site).real();
        result[1] *= Tr(1) - greenDs[split](site, site).real();
        result += Tr(1);
        return result; // Eq.(7.36) in [2]
    }

    template<Scalar T>
    auto DQMC<T>::calcLnSpinWaveWeight(int site) const noexcept -> Tr {
        assert(0 <= site && site < getNumSite());
        const auto spins = aux.row(site);
        return calcLnSpinWaveWeight(spins.sum());
    }

    template<Scalar T>
    auto DQMC<T>::calcLnSpinWaveWeight(Tr sumSpin) const noexcept -> Tr {
        Tr split = getNumSplit();
        Tr a = (split + sumSpin) * 0.5;
        Tr b = (split - sumSpin) * 0.5;
        Tr combine;
        if (a.isPositive() && b.isPositive())
            combine = (a + 0.5) * ln(a) + (b + 0.5) * ln(b);
        else
            combine = (split + 0.5) * ln(split);
        return ln1pexp(lncosh(params.getAlpha() * sumSpin) - params.getLnCoshShift()) - combine;
    }

    template<Scalar T>
    template<RNG R>
    bool DQMC<T>::accept(Tr deltaW) noexcept {
        if (deltaW.isPositive())
            return true;

        const auto p = Tr::template random_uniform<R>();
        return p < exp(deltaW);
    }

    template<Scalar T>
    bool DQMC<T>::accept(Tr deltaW, Tr prop) noexcept {
        if (deltaW.isPositive())
            return true;
        return prop < exp(deltaW);
    }
}
