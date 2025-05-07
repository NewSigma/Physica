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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/QRDecomp.h"
#include "Physica/Core/Physics/ManyBody/Model/Hubbard.h"

namespace Physica {
    template<Scalar T, int Dim, int NumSite = Dynamic>
    class DQMC {
        static_assert(!Diffable<T>);
        static_assert(!T::isComplex, "[Error]: Model param must be real");
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        using This = DQMC<T, Dim, NumSite>;
        using Tc = T::ComplexType;
        using Tf = Diff<T, DiffMode::Forward, 1>;
        using ModelType = Hubbard<T, Dim>;
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, NumSite, NumSite>;
    private:
        DenseSymmMatrix<T, NumSite> hoppingMatrix;
        ModelType hubbard;
        T alpha;
        T beta;
        T chemMu;

        VectorND<T> aux;

        QRDecomp<T> qr;
        Array<MatrixType> chain;
        MatrixType expB;
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
        DQMC() = default;
        DQMC(ModelType hubbard_, T beta_, T chemMu_, int numSplit);
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
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return hoppingMatrix; }
        [[nodiscard]] size_t getNumSite() const noexcept { return hubbard.getNumSuperCellSite(); }
        [[nodiscard]] T getHoppingT() const noexcept { return hubbard.getHoppingT(); }
        [[nodiscard]] T getRepelU() const noexcept { return hubbard.getRepelU(); }
        [[nodiscard]] T getBeta() const noexcept { return beta; }
        [[nodiscard]] const auto getAuxField() const noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }
        [[nodiscard]] int getNumSplit() const noexcept { return chain.getLength(); }

        [[nodiscard]] T getLnPartitionZ() const noexcept { return lnPartitionZ; }
        [[nodiscard]] T getSignU() const noexcept { return signU; }
        [[nodiscard]] T getSignD() const noexcept { return signD; }
        [[nodiscard]] T getSign() const noexcept { return signU * signD; }
        [[nodiscard]] const auto& getGreenU() const noexcept { return greenU; }
        [[nodiscard]] const auto& getGreenD() const noexcept { return greenD; }
        /* Setters */
        void setBeta(T beta_);
        void setChemMu(T chemMu_);
    private:
        void makeHoppingMatrix();
        void makeWeightMatrix(bool spin);
        std::pair<T, T> lnPartition(bool spin);
        void calcGreen(bool spin);

        T calcLnSpinWaveWeight(int split) const;
        T calcLnSpinWaveWeight(T sumSpin) const;
        /* Getters */
        [[nodiscard]] auto getAuxField() noexcept { return aux.reshape_col(getNumSite(), getNumSplit()); }
    };

    template<Scalar T, int Dim, int NumSite>
    DQMC<T, Dim, NumSite>::DQMC(ModelType hubbard_, T beta_, T chemMu_, int numSplit)
            : hubbard(std::move(hubbard_))
            , chain(numSplit) {
        assert(!hubbard.getRepelU().isNegative() && "[Error]: It is assumed U >= 0");
        assert(numSplit > 0 && "[Error]: Invalid NumSplit");
        const int numSite = getNumSite();
        qr.resize(numSite, numSite);
        aux.resize(numSite * numSplit);
        matrixDb.resize(numSite);
        matrixDs.resize(numSite);
        diagB1.resize(numSite);
        diagS1.resize(numSite);
        for (auto& m : chain)
            m.resize(numSite);
        greenU.resize(numSite, numSite);
        greenD.resize(numSite, numSite);

        makeHoppingMatrix();
        setBeta(std::move(beta_));
        setChemMu(std::move(chemMu_));
    }

    template<Scalar T, int Dim, int NumSite>
    template<RNG R>
    void DQMC<T, Dim, NumSite>::step() {
        const Array<int> sites = R::random_int(getNumSplit(), 0, getNumSite() - 1);
        auto field = getAuxField();
        for (int split = 0; split < getNumSplit(); ++split) {
            auto spins = field.col(split);
            const int site = sites[split];
            const T sumSpin0 = spins.sum();
            const T sumSpin1 = sumSpin0 - T(2) * spins[site];
            const T delta = calcLnSpinWaveWeight(sumSpin1) - calcLnSpinWaveWeight(sumSpin0);
            if (delta.isNegative()) {
                const T p = T::template random_uniform<R>();
                const bool accept = p < exp(delta);
                if (!accept)
                    return;
            }
            spins[site] = -spins[site];
        }

        calcGreen(true);
        calcGreen(false);
    }

    template<Scalar T, int Dim, int NumSite>
    template<RNG R>
    void DQMC<T, Dim, NumSite>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int i = 0; i < numStep; ++i)
            step<R>();
    }

    template<Scalar T, int Dim, int NumSite>
    template<RNG R>
    void DQMC<T, Dim, NumSite>::random_uniform() {
        aux.template random_uniform<R>();
        aux = unit(aux - T(0.5));
    }

    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hoppingMatrix.swap(obj.hoppingMatrix);
        hubbard.swap(obj.hubbard);
        alpha.swap(obj.alpha);
        beta.swap(obj.beta);
        chemMu.swap(obj.chemMu);

        aux.swap(obj.aux);

        qr.swap(obj.qr);
        chain.swap(obj.chain);
        expB.swap(obj.expB);
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

    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::setBeta(T beta_) {
        assert(!beta_.isNegative() && "[Error]: Negative temperature is invalid");
        beta = std::move(beta_);

        const T betaM = beta / T(getNumSplit());
        const T x = betaM * getRepelU();
        alpha = x * T(0.5) + ln1p(sqrt(T(1) - exp(-x)));

        DenseSymmMatrix<T, NumSite> hoppingMatrixB = hoppingMatrix * betaM;
        expB = exp(hoppingMatrixB);
    }

    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::setChemMu(T chemMu_) {
        chemMu = std::move(chemMu_);
    }

    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::makeHoppingMatrix() {
        const size_t numSite = getNumSite();
        hoppingMatrix.resize(numSite);
        hoppingMatrix = T(0);

        for (size_t from = 0; from < numSite; ++from) {
            if constexpr (ModelType::UntrivialNearestNeighbor) {
                const auto& targets = hubbard.getHopIndexArray()[from];
                for (size_t to : targets)
                    hoppingMatrix(from, to) = -getHoppingT();
            }
            else
                hoppingMatrix(from, (from + 1) % numSite) = -getHoppingT();
        }
    }
    /**
     * Reference:
     * [1] Linear Algebra and its Applications 435(3), 659-673 (2011); https://doi.org/10.1016/j.laa.2010.06.023
     */
    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::makeWeightMatrix(bool spin) {
        const int numSplit = getNumSplit();
        /* Init chain */ {
            const auto field = getAuxField();
            for (int i = 0; i < numSplit; ++i) {
                const auto col = field.col(i);
                auto& m = chain[i];
                m = expB;
                for (size_t j = 0; j < getNumSite(); ++j) {
                    const auto sigma = spin ? col.calc(j) : -col.calc(j);
                    m.col(j) *= exp(alpha * sigma);
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

        const T shift = beta * (chemMu - getRepelU() * T(0.5));
        for (size_t i = 0; i < getNumSite(); ++i) {
            const T lnD = diagB[i] + shift;
            const T expD = exp(lnD);
            const T sign = diagS[i];
            const bool sep = lnD.isPositive();
            diagB[i] = sep ? expD : T(1);
            diagS[i] = (sep ? T(1) : std::max(expD, T(std::numeric_limits<T>::min()))) * sign;
        }
    }

    template<Scalar T, int Dim, int NumSite>
    auto DQMC<T, Dim, NumSite>::lnPartition(bool spin) -> std::pair<T, T> {
        makeWeightMatrix(spin);
        T sign = qr.calcDetQ();
        buffer = qr.getMatrixQ() * matrixDb.inverse();
        qr.compute(buffer.transpose() + matrixDs * matrixR);

        T lnZ = matrixDb.lnAbsDet() + qr.getMatrixR().lnAbsDet();
        for (int i = 0; i < getNumSplit(); ++i)
            lnZ -= calcLnSpinWaveWeight(i);
        sign *= qr.calcDetQ() * unit(qr.getMatrixR().diag()).prod();
        return std::make_pair(lnZ, sign);
    }

    template<Scalar T, int Dim, int NumSite>
    void DQMC<T, Dim, NumSite>::calcGreen(bool spin) {
        const auto [lnZ, sign] = lnPartition(spin);
        lnPartitionZ = lnZ;
        if (spin)
            signU = sign;
        else
            signD = sign;

        auto& green = spin ? greenU : greenD;
        matrixQ = buffer * qr.getMatrixQ();
        green = qr.getMatrixR().inverse() * matrixQ.transpose();
    }

    template<Scalar T, int Dim, int NumSite>
    T DQMC<T, Dim, NumSite>::calcLnSpinWaveWeight(int split) const {
        assert(0 <= split && split < getNumSplit());
        const auto field = getAuxField();
        const auto spins = field.col(split);
        return calcLnSpinWaveWeight(spins.sum());
    }

    template<Scalar T, int Dim, int NumSite>
    T DQMC<T, Dim, NumSite>::calcLnSpinWaveWeight(T sumSpin) const {
        const T factor = lncosh(beta * (chemMu - getRepelU() * 0.5));
        return ln1pexp(lncosh(alpha * sumSpin) - factor);
    }
}
