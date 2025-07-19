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
     */
    template<Scalar T>
    class DQMC {
        using This = DQMC<T>;
        using Params = HubbardParams<T>;
        using MatrixND = Params::MatrixND;
    private:
        const Params& params;
        MatrixND aux;
        CyclicChainQDT<T> chainU;
        CyclicChainQDT<T> chainD;

        QRDecomp<T> qr;
        DiagMatrix<T> diagB;
        DiagMatrix<T> diagS;
        MatrixND buffer;
        T lnZ0;

        T lnPartitionZ;
        T signU = 1;
        T signD = 1;
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

        void update();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params.getNumSplit(); }
        [[nodiscard]] auto& getAuxField() noexcept { return aux; }

        [[nodiscard]] T getLnPartitionZ() const noexcept { return lnPartitionZ; }
        [[nodiscard]] T getSignU() const noexcept { return signU; }
        [[nodiscard]] T getSignD() const noexcept { return signD; }
        [[nodiscard]] T getSign() const noexcept { return signU * signD; }
        [[nodiscard]] int getNumEqualGreen() const noexcept { return greenUs.getLength(); }
        [[nodiscard]] int getPeriod() const noexcept { return getNumSplit() / getNumEqualGreen(); }
        [[nodiscard]] const auto& getGreenUs() const noexcept { return greenUs; }
        [[nodiscard]] const auto& getGreenDs() const noexcept { return greenDs; }
    private:
        /* Operations */
        void resize(int numSite, int numSplit, int period);

        void initChain();
        void single_flip(int site, int split) noexcept;
        const QDTDecomp<T>& makeWeightMatrix(CyclicChainQDT<T>& chain, int shift);
        std::pair<T, T> lnPartition(CyclicChainQDT<T>& chain, int shift);
        std::pair<T, T> lnPartition(int shift);
        std::pair<T, T> calcGreen(CyclicChainQDT<T>& chain, MatrixND& green, int shift);
        void calcGreens();

        T calcLnSpinWaveWeight(int site) const noexcept;
        T calcLnSpinWaveWeight(T sumSpin) const noexcept;
        /* Static members */
        template<RNG R>
        static bool accept(T deltaW) noexcept;
        static bool accept(T deltaW, T prop) noexcept;
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_, int period) : params(params_) {
        resize(params.getNumSite(), params.getNumSplit(), period);
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_random() {
        aux.template random_uniform<R>();
        aux = unit_elem(aux - T(0.5));
        initChain();
        calcGreens();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_mh() {
        auto& rng = R::getInstance();
        const int site = std::uniform_int_distribution<>(0, getNumSite() - 1)(rng);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(rng);
        single_flip(site, split);

        const T lnZ1 = lnPartition().first;
        if (!accept<R>(lnZ1 - lnZ0)) {
            single_flip(site, split);
            return;
        }

        calcGreens();
        lnZ0 = lnPartitionZ;
        lnPartitionZ = 0;
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_mh_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        auto props = VectorND<T>(getNumSite());
        for (int step = 0, split = 0; step < numStep; ++step) {
            props.template random_uniform<R>();
            for (int site = 0; site < getNumSite(); ++site) {
                single_flip(site, split);

                const T lnZ1 = lnPartition().first;
                if (!accept(lnZ1 - lnZ0, props[site])) {
                    single_flip(site, split);
                    continue;
                }
                lnPartitionZ = lnZ1;
            }
            split = (split + 1) % getNumSplit();
        }
        calcGreens();
        lnZ0 = lnPartitionZ;
        lnPartitionZ = 0;
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin() {
        auto& rng = R::getInstance();
        const int site = std::uniform_int_distribution<>(0, getNumSite() - 1)(rng);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(rng);
        auto spins = aux.row(site);
        const T sumSpin0 = spins.sum();
        const T sumSpin1 = sumSpin0 - T(2) * spins[split];
        const T deltaW = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
        const bool flip1 = accept<R>(deltaW);
        if (flip1)
            single_flip(site, split);

        const T p = T::template random_uniform<R>();
        const bool flip2 = p < 0.5;
        if (flip2) {
            for (int i = 0; i < getNumSplit(); ++i)
                single_flip(site, i);
        }

        const bool noFlip = !flip1 && !flip2;
        if (noFlip)
            return;
    
        calcGreens();
        for (int i = 0; i < getNumSite(); ++i)
            lnPartitionZ -= calcLnSpinWaveWeight(i);
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_spin_for(int numStep)  {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        assert(numStep % getNumSplit() == 0 && "[Warn]: Suggest divisable step num");
        auto props = VectorND<T>(getNumSplit());
        auto splits = Array<int>(getNumSplit());
        for (int i = 0; i < getNumSplit(); ++i)
            splits[i] = i;

        for (int site = 0; site < getNumSite(); ++site) {
            auto spins = aux.row(site);
            T sumSpin0 = spins.sum();
            T weight0 = calcLnSpinWaveWeight(sumSpin0);
            for (int step = 0; step < numStep; ++step) {
                const int i = step % getNumSplit();
                if (i == 0) {
                    props.template random_uniform<R>();
                    std::shuffle(splits.begin(), splits.end(), R::getInstance());
                }

                const int split = splits[i];
                const T sumSpin1 = sumSpin0 - T(2) * spins[split];
                const T weight1 = calcLnSpinWaveWeight(sumSpin1);
                const T deltaW = weight1 - weight0;
                if (accept(deltaW, props[i])) {
                    spins[split] = -spins[split];
                    sumSpin0 = sumSpin1;
                    weight0 = weight1;
                }
            }

            const T p = T::template random_uniform<R>();
            if (p < 0.5)
                spins = -spins;
        }
        initChain();
        calcGreens();
        for (int i = 0; i < getNumSite(); ++i)
            lnPartitionZ -= calcLnSpinWaveWeight(i);
    }

    template<Scalar T>
    void DQMC<T>::update() {
        initChain();
        calcGreens();
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        aux.swap(obj.aux);

        chainU.swap(obj.chainU);
        chainD.swap(obj.chainD);

        qr.swap(obj.qr);
        diagB.swap(obj.diagB);
        diagS.swap(obj.diagS);
        buffer.swap(obj.buffer);

        lnZ0.swap(obj.lnZ0);
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

        qr.resize(numSite, numSite);
        diagB.resize(numSite);
        diagS.resize(numSite);
        buffer.resize(numSite, numSite);

        int numEqualGreen = numSplit / period;
        greenUs.resize(numEqualGreen, numSite, numSite);
        greenDs.resize(numEqualGreen, numSite, numSite);
    }

    template<Scalar T>
    void DQMC<T>::initChain() {
        const int numSplit = getNumSplit();
        const auto& expB = params.getExpB();
        DiagMatrix<T> expU(getNumSite());
        for (int split = 0; split < numSplit; ++split) {
            expU.diag() = exp(params.getAlpha() * aux.col(split));
            chainU[split].compute(expB * expU);
            expU.diag() = exp(-params.getAlpha() * aux.col(split));
            chainD[split].compute(expB * expU);
        }
        chainU.invalidates();
        chainD.invalidates();
    }

    template<Scalar T>
    void DQMC<T>::single_flip(int site, int split) noexcept {
        auto spins = aux.row(site);
        spins[split] = -spins[split];

        DiagMatrix<T> expU(getNumSite());
        const auto& expB = params.getExpB();
        expU.diag() = exp(params.getAlpha() * aux.col(split));
        chainU[split].compute(expB * expU);
        expU.diag() = exp(-params.getAlpha() * aux.col(split));
        chainD[split].compute(expB * expU);
        chainU.invalidate(split);
        chainD.invalidate(split);
    }

    template<Scalar T>
    const QDTDecomp<T>& DQMC<T>::makeWeightMatrix(CyclicChainQDT<T>& chain, int shift) {
        const int numSplit = getNumSplit();
        const auto& decomp = chain.multiply(shift, (numSplit + shift - 1) % numSplit);
        const auto& matrixD = decomp.getMatrixD();
        const T betaMu = params.calcBetaMu();
        for (int i = 0; i < getNumSite(); ++i) {
            const T originD = matrixD.diag()[i] * exp(betaMu);
            const T absD = abs(originD);
            const bool sep = absD > T(1);
            diagB.diag()[i] = sep ? reciprocal(absD) : T(1);
            diagS.diag()[i] = sep ? unit(originD) : originD;
        }
        return decomp;
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition(CyclicChainQDT<T>& chain, int shift) {
        const auto& qdt = makeWeightMatrix(chain, shift);
        buffer = qdt.getMatrixQ() * diagB;
        qr.compute(buffer.transpose() + diagS * qdt.getMatrixT());
        // Handle potential underflow
        diagB.diag() += T(std::numeric_limits<T>::min());
        qr.getWorking().diag() += T(std::numeric_limits<T>::min());

        T lnZ = -diagB.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        T sign = qdt.calcDetQ() * qr.calcDetQ() * unit(qr.getMatrixR().diag()).prod();
        return std::make_pair(lnZ, sign);
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition(int shift) {
        auto [lnZ1, sign1] = lnPartition(chainU, shift);
        auto [lnZ2, sign2] = lnPartition(chainD, shift);
        return std::make_pair(lnZ1 + lnZ2, sign1 * sign2);
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::calcGreen(CyclicChainQDT<T>& chain, MatrixND& green, int shift) {
        const auto pair = lnPartition(chain, shift);
        MatrixND temp = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inverse() * temp.transpose();
        return pair;
    }

    template<Scalar T>
    void DQMC<T>::calcGreens() {
        const int period = getPeriod();
        for (int shift = 0, i = 0; shift < getNumSplit(); shift += period, ++i) {
            auto [lnZ1, sign1] = calcGreen(chainU, greenUs[i], shift);
            assert((shift == 0 || signU == sign1) && "[Error]: Unexpected sign mismatch");
            signU = sign1;

            auto [lnZ2, sign2] = calcGreen(chainD, greenDs[i], shift);
            assert((shift == 0 || signD == sign2) && "[Error]: Unexpected sign mismatch");
            signD = sign2;

            toNextMean(lnPartitionZ, i, lnZ1 + lnZ2);
        }
    }

    template<Scalar T>
    T DQMC<T>::calcLnSpinWaveWeight(int site) const noexcept {
        assert(0 <= site && site < getNumSite());
        const auto spins = aux.row(site);
        return calcLnSpinWaveWeight(spins.sum());
    }

    template<Scalar T>
    T DQMC<T>::calcLnSpinWaveWeight(T sumSpin) const noexcept {
        T split = getNumSplit();
        T a = (split + sumSpin) * 0.5;
        T b = (split - sumSpin) * 0.5;
        T combine;
        if (a.isPositive() && b.isPositive())
            combine = (a + 0.5) * ln(a) + (b + 0.5) * ln(b);
        else
            combine = (split + 0.5) * ln(split);
        return ln1pexp(lncosh(params.getAlpha() * sumSpin) - params.getLnCoshShift()) - combine;
    }

    template<Scalar T>
    template<RNG R>
    bool DQMC<T>::accept(T deltaW) noexcept {
        if (deltaW.isPositive())
            return true;

        const T p = T::template random_uniform<R>();
        return p < exp(deltaW);
    }

    template<Scalar T>
    bool DQMC<T>::accept(T deltaW, T prop) noexcept {
        if (deltaW.isPositive())
            return true;
        return prop < exp(deltaW);
    }
}
