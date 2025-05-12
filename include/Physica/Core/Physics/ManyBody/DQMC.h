/*
 * Copyright 2025 WeiBo He.
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
    template<Scalar T>
    class DQMC {
        using This = DQMC<T>;
        using Params = HubbardParams<T>;
    public:
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
    private:
        const Params& params;
        VectorND<T> aux;
        QRDecomp<T> qr;
        Array<MatrixType> chain;

        DiagMatrix<T> matrixDb;
        DiagMatrix<T> matrixDs;
        VectorND<T> diagB1;
        VectorND<T> diagS1;
        MatrixType matrixQ;
        MatrixType matrixR;
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
        template<RNG R>
        void step();
        template<RNG R>
        void step_for(int numStep);

        template<RNG R>
        void random_uniform();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return qr.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return chain.getLength(); }
        [[nodiscard]] const auto getAuxField() const noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }
        [[nodiscard]] auto getAuxField() noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }

        [[nodiscard]] T getLnPartitionZ() const noexcept { return lnPartitionZ; }
        [[nodiscard]] T getSignU() const noexcept { return signU; }
        [[nodiscard]] T getSignD() const noexcept { return signD; }
        [[nodiscard]] T getSign() const noexcept { return signU * signD; }
        [[nodiscard]] const auto& getGreenU() const noexcept { return greenU; }
        [[nodiscard]] const auto& getGreenD() const noexcept { return greenD; }
    private:
        DQMC(const Params& params_, int numSite, int numSplit);
        void makeWeightMatrix(bool spin);
        std::pair<T, T> lnPartition(bool spin);
        T calcGreen(bool spin);
        void update();

        T calcLnSpinWaveWeight(int site) const;
        T calcLnSpinWaveWeight(T sumSpin) const;
    };

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_)
            : DQMC(params_, params_.getNumSite(), params_.getNumSplit()) {}

    template<Scalar T>
    DQMC<T>::DQMC(const Params& params_, int numSite, int numSplit)
            : params(params_)
            , aux(numSite * numSplit)
            , qr(numSite, numSite)
            , chain(numSplit)
            , matrixDb(numSite)
            , matrixDs(numSite)
            , diagB1(numSite)
            , diagS1(numSite)
            , greenU(numSite, numSite)
            , greenD(numSite, numSite) {
        assert(numSite > 0 && "[Error]: Invalid NumSite");
        assert(numSplit > 0 && "[Error]: Invalid NumSplit");
        for (auto& m : chain)
            m.resize(numSite);
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step() {
        const Array<int> splits = R::random_int(getNumSite(), 0, getNumSplit() - 1);
        const auto probs = VectorND<T>::template random_uniform<R>(getNumSplit());

        auto field = getAuxField();
        for (int site = 0; site < getNumSite(); ++site) {
            auto spins = field.row(site);
            spins *= (probs[site] > T(0.5)) ? T(1) : T(-1);

            const int split = splits[site];
            const T sumSpin0 = spins.sum();
            const T sumSpin1 = sumSpin0 - T(2) * spins[split];
            const T delta = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
            if (delta.isNegative()) {
                const T p = T::template random_uniform<R>();
                const bool accept = p < exp(delta);
                if (!accept)
                    continue;
            }
            spins[split] = -spins[split];
        }
        update();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int i = 0; i < numStep; ++i)
            step<R>();
    }

    template<Scalar T>
    template<RNG R>
    void DQMC<T>::random_uniform() {
        aux.template random_uniform<R>();
        aux = unit(aux - T(0.5));
        update();
    }

    template<Scalar T>
    void DQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        aux.swap(obj.aux);
        qr.swap(obj.qr);
        chain.swap(obj.chain);
        
        matrixDb.swap(obj.matrixDb);
        matrixDs.swap(obj.matrixDs);
        diagB1.swap(obj.diagB1);
        diagS1.swap(obj.diagS1);
        matrixQ.swap(obj.matrixQ);
        matrixR.swap(obj.matrixR);
        buffer.swap(obj.buffer);

        lnPartitionZ.swap(obj.lnPartitionZ);
        signU.swap(obj.signU);
        signD.swap(obj.signD);
        greenU.swap(obj.greenU);
        greenD.swap(obj.greenD);
    }
    /**
     * Reference:
     * [1] Linear Algebra and its Applications 435(3), 659-673 (2011); https://doi.org/10.1016/j.laa.2010.06.023
     */
    template<Scalar T>
    void DQMC<T>::makeWeightMatrix(bool spin) {
        const int numSplit = getNumSplit();
        /* Init chain */ {
            const auto field = getAuxField();
            for (int i = 0; i < numSplit; ++i) {
                const auto col = field.col(i);
                auto& m = chain[i];
                m = params.getExpB();
                for (size_t j = 0; j < getNumSite(); ++j) {
                    const auto sigma = spin ? col.calc(j) : -col.calc(j);
                    m.col(j) *= exp(params.getAlpha() * sigma);
                }
            }
        }

        auto& diagB = matrixDb.diag();
        auto& diagS = matrixDs.diag();
        qr.compute(chain[numSplit - 1]);

        diagB = qr.toQDT();
        diagS = unit(diagB);
        diagB = ln(abs(diagB));

        matrixR = qr.getMatrixR();
        for (int i = numSplit - 2; i >= 0; --i) {
            qr.compute(chain[i] * qr.getMatrixQ());
            /* Make new diag matrix */ {
                const auto triu = qr.getMatrixR();
                const auto diagR = triu.diag();
                diagB1 = diagB + ln(abs(diagR));
                diagS1 = unit(diagR);
            }
            auto& working = qr.getWorking();
            for (size_t r = 0; r < getNumSite(); ++r) {
                for (size_t c = r; c < getNumSite(); ++c) {
                    if (r == c) {
                        working(r, c) = 1;
                        continue;
                    }
                    working(r, c) *= -exp(-diagB1[r] + diagB[c]) * (diagS1[r] * diagS[c]);
                }
            }
            diagB1.swap(diagB);
            diagS = hadamard(diagS, diagS1);

            buffer = qr.getMatrixR() * matrixR;
            buffer.swap(matrixR);
        }

        const T shift = params.calcShift();
        for (size_t i = 0; i < getNumSite(); ++i) {
            const T lnD = diagB[i] + shift;
            const T expD = exp(lnD);
            const T sign = diagS[i];
            const bool sep = lnD.isPositive();
            diagB[i] = sep ? expD : T(1);
            diagS[i] = (sep ? T(1) : std::max(expD, T(std::numeric_limits<T>::min()))) * sign;
        }
    }

    template<Scalar T>
    auto DQMC<T>::lnPartition(bool spin) -> std::pair<T, T> {
        makeWeightMatrix(spin);
        T sign = qr.calcDetQ();
        buffer = qr.getMatrixQ() * matrixDb.inverse();
        qr.compute(buffer.transpose() + matrixDs * matrixR);

        T lnZ = matrixDb.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        sign *= qr.calcDetQ() * unit(qr.getMatrixR().diag()).prod();
        return std::make_pair(lnZ, sign);
    }

    template<Scalar T>
    T DQMC<T>::calcGreen(bool spin) {
        const auto [lnZ, sign] = lnPartition(spin);
        if (spin)
            signU = sign;
        else
            signD = sign;

        auto& green = spin ? greenU : greenD;
        matrixQ = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inverse() * matrixQ.transpose();
        return lnZ;
    }

    template<Scalar T>
    void DQMC<T>::update() {
        lnPartitionZ = calcGreen(true) + calcGreen(false);
        for (int i = 0; i < getNumSite(); ++i)
            lnPartitionZ -= calcLnSpinWaveWeight(i);
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
}
