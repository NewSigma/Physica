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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/QRDecomp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "DQMCImpl/HubbardParams.h"

namespace Physica {
    /**
     * Reference:
     * [1] Phys. Rev. B 40, 506; https://doi.org/10.1103/PhysRevB.40.506
     */
    template<Scalar T>
    class DQMC {
        using This = DQMC<T>;
        using Params = HubbardParams<T>;
        using MatrixType = Params::MatrixType;
        using MatrixChain = Array<MatrixType>;
    public:
        enum FlipMethod {
            MH,
            Spin,
            Uniform
        };
    private:
        const Params& params;
        VectorND<T> aux;
        MatrixChain chainU;
        MatrixChain chainD;

        Array<QRDecomp<T>> qr;
        Array<MatrixType> matrixQs;
        Array<DiagMatrix<T>> matrixDs;
        MatrixType buffer;
        T lnZ0;

        T lnPartitionZ;
        T signU = 1;
        T signD = 1;
        MatrixType greenU;
        MatrixType greenD;
    public:
        DQMC() = delete;
        DQMC(const Params& params_);
        DQMC(const This&) = default;
        DQMC(This&&) noexcept = default;
        ~DQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R, FlipMethod Method>
        void step();
        template<RNG R, FlipMethod Method>
        void step_for(int numStep);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params.getNumSplit(); }
        [[nodiscard]] const auto getAuxField() const noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }
        [[nodiscard]] auto getAuxField() noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }

        [[nodiscard]] T getLnPartitionZ() const noexcept { return lnPartitionZ; }
        [[nodiscard]] T getSignU() const noexcept { return signU; }
        [[nodiscard]] T getSignD() const noexcept { return signD; }
        [[nodiscard]] T getSign() const noexcept { return signU * signD; }
        [[nodiscard]] const auto& getGreenU() const noexcept { return greenU; }
        [[nodiscard]] const auto& getGreenD() const noexcept { return greenD; }
    private:
        /* Operations */
        void resize(int numSite, int numSplit);
        template<RNG R>
        void random_uniform();

        void initChain();
        void single_flip(int site, int split);
        void makeWeightMatrix(const MatrixChain& chain);
        std::pair<T, T> lnPartition(const MatrixChain& chain);
        std::pair<T, T> lnPartition();
        std::pair<T, T> calcGreen(const MatrixChain& chain, MatrixType& green);
        template<FlipMethod Method>
        void update();

        T calcLnSpinWaveWeight(int site) const;
        T calcLnSpinWaveWeight(T sumSpin) const;
        /* Getters */
        [[nodiscard]] auto& getDiagB() noexcept { return matrixDs[0]; }
        [[nodiscard]] auto& getDiagS() noexcept { return matrixDs[1]; }
        /* Static members */
        template<RNG R>
        static bool accept(T deltaW) noexcept;
        static bool accept(T deltaW, T prop) noexcept;
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_) : params(params_) {
        resize(params.getNumSite(), params.getNumSplit());
    }

    template<Scalar T>
    template<RNG R, DQMC<T>::FlipMethod Method>
    void DQMC<T>::step() {
        auto& rng = R::getInstance();
        const int site = std::uniform_int_distribution<>(0, getNumSite() - 1)(rng);
        const int split = std::uniform_int_distribution<>(0, getNumSplit() - 1)(rng);
        if constexpr (Method == MH) {
            single_flip(site, split);

            const T lnZ1 = lnPartition().first;
            if (!accept<R>(lnZ1 - lnZ0)) {
                single_flip(site, split);
                return;
            }
        }
        else if constexpr (Method == Spin) {
            auto field = getAuxField();
            auto spins = field.row(site);
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
        }
        else
            random_uniform<R>();
        update<Method>();
    }

    template<Scalar T>
    template<RNG R, DQMC<T>::FlipMethod Method>
    void DQMC<T>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        if constexpr (Method == Uniform) {
            random_uniform<R>();
            update<Method>();
            return;
        }

        auto splits = Array<int>(getNumSite());
        auto props = VectorND<T>(getNumSite());
        if constexpr (Method == MH) {
            for (int step = 0; step < numStep; ++step) {
                R::random_int(splits, 0, getNumSplit() - 1);
                props.template random_uniform<R>();
                for (int site = 0; site < getNumSite(); ++site) {
                    const int split = splits[site];
                    single_flip(site, split);

                    const T lnZ1 = lnPartition().first;
                    if (!accept(lnZ1 - lnZ0, props[site])) {
                        single_flip(site, split);
                        continue;
                    }
                    lnPartitionZ = lnZ1;
                }
            }
        }
        else if constexpr (Method == Spin) {
            for (int step = 0; step < numStep; ++step) {
                R::random_int(splits, 0, getNumSplit() - 1);
                props.template random_uniform<R>();
                for (int site = 0; site < getNumSite(); ++site) {
                    const int split = splits[site];
                    auto field = getAuxField();
                    auto spins = field.row(site);
                    const T sumSpin0 = spins.sum();
                    const T sumSpin1 = sumSpin0 - T(2) * spins[split];
                    const T deltaW = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
                    if (accept(deltaW, props[site]))
                        spins[split] = -spins[split];
                }
            }

            for (int site = 0; site < getNumSite(); ++site) {
                const T p = T::template random_uniform<R>();
                auto field = getAuxField();
                auto spins = field.row(site);
                if (p < 0.5)
                    spins = -spins;
            }
            initChain();
        }
        else
            noImpl("Unknown Flip Method");
        update<Method>();
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        aux.swap(obj.aux);

        chainU.swap(obj.chainU);
        chainD.swap(obj.chainD);

        qr.swap(obj.qr);
        matrixQs.swap(obj.matrixQs);
        matrixDs.swap(obj.matrixDs);
        buffer.swap(obj.buffer);

        lnPartitionZ.swap(obj.lnPartitionZ);
        signU.swap(obj.signU);
        signD.swap(obj.signD);
        greenU.swap(obj.greenU);
        greenD.swap(obj.greenD);
    }

    template<Scalar T>
    void DQMC<T>::resize(int numSite, int numSplit) {
        assert(numSite > 0 && "[Error]: Invalid NumSite");
        assert(numSplit > 0 && "[Error]: Invalid NumSplit");
        const int numQR = (numSplit + 1) / 2;
        aux.resize(numSite * numSplit);
        chainU.resize(numSplit);
        chainD.resize(numSplit);
        for (int i = 0; i < numSplit; ++i) {
            chainU[i].resize(numSite);
            chainD[i].resize(numSite);
        }

        qr.resize(numQR, numSite, numSite);
        matrixQs.resize(numQR, numSite, numSite);
        matrixDs.resize(numQR, numSite);
        buffer.resize(numSite, numSite);
        greenU.resize(numSite, numSite);
        greenD.resize(numSite, numSite);
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::random_uniform() {
        aux.template random_uniform<R>();
        aux = unit(aux - T(0.5));
        initChain();
    }

    template<Scalar T>
    void DQMC<T>::initChain() {
        const int numSplit = getNumSplit();
        const auto field = getAuxField();
        const auto& expB = params.getExpB();
        for (int split = 0; split < numSplit; ++split) {
            auto& mU = chainU[split];
            auto& mD = chainD[split];
            for (int site = 0; site < getNumSite(); ++site) {
                const T factor = params.getAlpha() * field(site, split);
                mU.col(site) = expB.col(site) * exp(factor);
                mD.col(site) = expB.col(site) * exp(-factor);
            }
        }
    }

    template<Scalar T>
    void DQMC<T>::single_flip(int site, int split) {
        auto field = getAuxField();
        auto spins = field.row(site);
        spins[split] = -spins[split];

        const T factor = T(2) * params.getAlpha() * spins[split];
        chainU[split].col(site) *= exp(factor);
        chainD[split].col(site) *= exp(-factor);
    }
    /**
     * Reference:
     * [1] Linear Algebra and its Applications 435(3), 659-673 (2011); https://doi.org/10.1016/j.laa.2010.06.023
     */
    template<Scalar T>
    void DQMC<T>::makeWeightMatrix(const MatrixChain& chain) {
        const int numSplit = getNumSplit();
        for (int i = 0; i < numSplit; i += 2) {
            const int index = i / 2;
            auto& qri = qr[index];
            if (i + 1 < numSplit) [[likely]]
                qri.compute(chain[i] * chain[i + 1]);
            else
                qri.compute(chain[i]);
            matrixQs[index] = qri.getMatrixQ();
            qri.toQDT(matrixDs[index].diag());
            qri.getWorking() = qri.getMatrixR();
        }

        int step = 1;
        while (step < qr.getLength()) {
            for (int i = 0; i + step < qr.getLength(); i += step * 2) {
                const int i1 = i;
                const int i2 = i + step;
                auto& d1 = matrixDs[i1];
                auto& d2 = matrixDs[i2];

                auto& qr1 = qr[i1];
                auto& qr2 = qr[i2];
                buffer = qr1.getMatrixR() * matrixQs[i2];
                buffer = d1 * buffer;
                buffer = buffer * d2;

                qr1.compute(buffer);
                qr1.toQDT(d1.diag());
                buffer = matrixQs[i1] * qr1.getMatrixQ();
                buffer.swap(matrixQs[i1]);
                buffer = qr1.getMatrixR() * qr2.getMatrixR();
                buffer.swap(qr1.getWorking());
            }
            step *= 2;
        }

        const T shift = params.calcShift();
        const auto& matrixD0 = matrixDs[0];
        for (int i = 0; i < getNumSite(); ++i) {
            const T originD = matrixD0.diag()[i] * exp(shift);
            const T absD = abs(originD);
            const bool sep = absD > T(1);
            getDiagB().diag()[i] = sep ? reciprocal(absD) : T(1);
            getDiagS().diag()[i] = sep ? originD.unit() : originD;
        }
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition(const MatrixChain& chain) {
        makeWeightMatrix(chain);
        auto& qr0 = qr[0];
        T sign = qr0.calcDetQ();
        buffer = matrixQs[0] * getDiagB();
        qr0.compute(buffer.transpose() + getDiagS() * qr0.getWorking());
        // Handle potential underflow
        getDiagB().diag() += T(std::numeric_limits<T>::min());
        qr0.getWorking().diag() += T(std::numeric_limits<T>::min());

        T lnZ = -getDiagB().lnAbsDet() + qr0.getMatrixR().lnAbsDet();
        sign *= qr0.calcDetQ() * unit(qr0.getMatrixR().diag()).prod();
        return std::make_pair(lnZ, sign);
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition() {
        auto [lnZ1, sign1] = lnPartition(chainU);
        auto [lnZ2, sign2] = lnPartition(chainD);
        return std::make_pair(lnZ1 + lnZ2, sign1 * sign2);
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::calcGreen(const MatrixChain& chain, MatrixType& green) {
        const auto pair = lnPartition(chain);
        const auto& qr0 = qr[0];
        matrixQs[0] = buffer * qr0.getMatrixQ();
        green = qr0.getMatrixR().inverse() * matrixQs[0].transpose();
        return pair;
    }

    template<Scalar T>
    template<DQMC<T>::FlipMethod Method>
    void DQMC<T>::update() {
        auto [lnZ1, sign1] = calcGreen(chainU, greenU);
        auto [lnZ2, sign2] = calcGreen(chainD, greenD);
        lnPartitionZ = lnZ1 + lnZ2;
        if constexpr (Method == MH) {
            lnZ0 = lnPartitionZ;
            lnPartitionZ = 0;
        }
        else if constexpr (Method == Spin) {
            for (int i = 0; i < getNumSite(); ++i)
                lnPartitionZ -= calcLnSpinWaveWeight(i);
        }
        signU = sign1;
        signD = sign2;
    }

    template<Scalar T>
    T DQMC<T>::calcLnSpinWaveWeight(int site) const {
        assert(0 <= site && site < getNumSite());
        const auto field = getAuxField();
        const auto spins = field.row(site);
        return calcLnSpinWaveWeight(spins.sum());
    }

    template<Scalar T>
    T DQMC<T>::calcLnSpinWaveWeight(T sumSpin) const {
        return ln1pexp(lncosh(params.getAlpha() * sumSpin) - params.getLnCoshShift());
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
