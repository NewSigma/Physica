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
        QRDecomp<T> qr;
        MatrixChain chainU;
        MatrixChain chainD;

        DiagMatrix<T> matrixDb;
        DiagMatrix<T> matrixDs;
        MatrixType matrixT;
        MatrixType buffer;

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
        void makeWeightMatrix(bool spin);
        std::pair<T, T> lnPartition(bool spin);
        std::pair<T, T> lnPartition();
        T calcGreen(bool spin);
        void update();

        T calcLnSpinWaveWeight(int site) const;
        T calcLnSpinWaveWeight(T sumSpin) const;
        /* Static members */
        template<RNG R>
        static bool accept(T deltaW) noexcept;
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

            const T lnZ0 = lnPartitionZ;
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
                for (int i = 0; i < getNumSite(); ++i)
                    single_flip(i, split);
            }

            const bool noFlip = !flip1 && !flip2;
            if (noFlip)
                return;
        }
        else
            random_uniform<R>();
        update();
    }

    template<Scalar T>
    template<RNG R, DQMC<T>::FlipMethod Method>
    void DQMC<T>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        Array<int> splits(getNumSite());
        if constexpr (Method == MH) {
            for (int step = 0; step < numStep; ++step) {
                R::random_int(splits, 0, getNumSplit() - 1);
                for (int site = 0; site < getNumSite(); ++site) {
                    const int split = splits[site];
                    single_flip(site, split);

                    const T lnZ0 = lnPartitionZ;
                    const T lnZ1 = lnPartition().first;
                    if (!accept<R>(lnZ1 - lnZ0)) {
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
                for (int site = 0; site < getNumSite(); ++site) {
                    auto field = getAuxField();
                    auto spins = field.row(site);
                    const int split = splits[site];
                    const T sumSpin0 = spins.sum();
                    const T sumSpin1 = sumSpin0 - T(2) * spins[split];
                    const T deltaW = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
                    if (accept<R>(deltaW))
                        spins[split] = -spins[split];

                    const T p = T::template random_uniform<R>();
                    if (p < 0.5)
                        spins = -spins;
                }
            }
            initChain();
        }
        else
            random_uniform<R>();
        update();
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        aux.swap(obj.aux);
        qr.swap(obj.qr);
        chainU.swap(obj.chainU);
        chainD.swap(obj.chainD);
        
        matrixDb.swap(obj.matrixDb);
        matrixDs.swap(obj.matrixDs);
        matrixT.swap(obj.matrixT);
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
        aux.resize(numSite * numSplit);
        qr.resize(numSite, numSite);
        chainU.resize(numSplit);
        chainD.resize(numSplit);
        for (int i = 0; i < numSplit; ++i) {
            chainU[i].resize(numSite);
            chainD[i].resize(numSite);
        }

        matrixDb.resize(numSite);
        matrixDs.resize(numSite);
        matrixT.resize(numSite, numSite);
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
    void DQMC<T>::makeWeightMatrix(bool spin) {
        const int numSplit = getNumSplit();
        const auto& chain = spin ? chainU : chainD;
        qr.compute(chain[numSplit - 1]);
        qr.toQDT();
        matrixT = qr.getMatrixR();
        for (int i = numSplit - 2; i >= 0; --i) {
            qr.compute(chain[i] * qr.getMatrixQ() * DiagMatrix<T>(qr.getVecD()));
            qr.toQDT();
            buffer = qr.getMatrixR() * matrixT;
            buffer.swap(matrixT);
        }

        const T shift = params.calcShift();
        for (int i = 0; i < getNumSite(); ++i) {
            const T originD = qr.getVecD()[i] * exp(shift);
            const T absD = abs(originD);
            const bool sep = absD > T(1);
            matrixDb.diag()[i] = sep ? reciprocal(absD) : T(1);
            matrixDs.diag()[i] = sep ? originD.unit() : originD;
        }
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition(bool spin) {
        makeWeightMatrix(spin);
        T sign = qr.calcDetQ();
        buffer = qr.getMatrixQ() * matrixDb;
        qr.compute(buffer.transpose() + matrixDs * matrixT);
        // Handle potential underflow
        matrixDb.diag() += T(std::numeric_limits<T>::min());
        qr.getWorking().diag() += T(std::numeric_limits<T>::min());

        T lnZ = -matrixDb.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        sign *= qr.calcDetQ() * unit(qr.getMatrixR().diag()).prod();
        return std::make_pair(lnZ, sign);
    }

    template<Scalar T>
    std::pair<T, T> DQMC<T>::lnPartition() {
        auto [lnZ1, sign1] = lnPartition(true);
        auto [lnZ2, sign2] = lnPartition(false);
        return std::make_pair(lnZ1 + lnZ2, sign1 * sign2);
    }

    template<Scalar T>
    T DQMC<T>::calcGreen(bool spin) {
        const auto [lnZ, sign] = lnPartition(spin);
        if (spin)
            signU = sign;
        else
            signD = sign;

        auto& green = spin ? greenU : greenD;
        matrixT = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inverse() * matrixT.transpose();
        return lnZ;
    }

    template<Scalar T>
    void DQMC<T>::update() {
        lnPartitionZ = calcGreen(true) + calcGreen(false);
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
        const T factor = lncosh(params.calcShift());
        return ln1pexp(lncosh(params.getAlpha() * sumSpin) - factor);
    }

    template<Scalar T>
    template<RNG R>
    bool DQMC<T>::accept(T deltaW) noexcept {
        if (deltaW.isPositive())
            return true;

        const T p = T::template random_uniform<R>();
        return p < exp(deltaW);
    }
}
