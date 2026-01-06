/*
 * Copyright 2025-2026 Weibo He.
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

#include "DQMCImpl/ActionMatrix.h"
#include "DQMCImpl/ImagKinetic.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"

namespace Physica {
    template<Scalar T>
    class FreqDQMC {
        using This = FreqDQMC<T>;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        using Cell = MDCell<Tr, 1>;
        using GreenPair = ImagKinetic<Tr>::GreenPair;
    private:
        Array<DenseLU<T, false>, 2> lu;
        ActionMatrix<T> action;

        GreenPair greens;
        Trv lnWeight = Trv::nan();
        Trv sign = 1;

        uint64_t numTotal = 0;
        uint64_t numAccept = 0;
    public:
        FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity);
        FreqDQMC(const This&) = default;
        FreqDQMC(This&&) noexcept = default;
        ~FreqDQMC() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R, ExecutePolicy P = Sequential>
        void step_random();
        template<RNG R, ExecutePolicy P = Sequential>
        void step(bool warmup = false);
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep);

        template<RNG R>
        void step_random(HamiltonMC<Tr>& hmc);
        template<RNG R, ExecutePolicy P = Sequential>
        void step(HamiltonMC<Tr>& hmc, bool warmup = false);
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep, HamiltonMC<Tr>& hmc);

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return action.getParams(); }
        [[nodiscard]] auto& getAuxField() noexcept { return action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return getParams().getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] const auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Trv getSign() const noexcept { return sign; }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        [[nodiscard]] Trv getAcceptRate() const noexcept { return Trv(numAccept) / Trv(numTotal); }
        /* Static members */
        [[nodiscard]] static int calcFreqCutoff(Trv beta, Trv freqDensity) noexcept;
    private:
        void makeFreqKernel(MatrixND<T>& freqKernel, const MatrixND<T>& invAction, int site) const;

        template<ExecutePolicy P>
        [[nodiscard]] Vector2D<Trv> calcDet();
        template<ExecutePolicy P>
        [[nodiscard]] Vector2D<Trv> calcLnWeight();
        template<ExecutePolicy P>
        auto calcGreen();
        /* Getters */
        [[nodiscard]] Trv getBetaU() const noexcept;
    };

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, calcFreqCutoff(params.getBeta(), freqDensity))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step_random() {
        action.template randAuxField<R>();
        auto [lnAD, sgnD] = calcLnWeight<P>();
        lnWeight = lnAD;
        sign = sgnD;
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step(bool warmup) {
        assert(lnWeight.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const int freq = std::uniform_int_distribution<int>(0, getNumFreq() * 2 - 1)(R::getInstance());
        const T save = action.template randAuxField<R>(freq, site);

        auto [lnW, sgnD] = calcLnWeight<P>();
        const bool accept = Trv::template random_uniform<R>() < exp(lnW - lnWeight);
        if (accept) {
            lnWeight = lnW;
            sign = sgnD;
            if (!warmup)
                calcGreen<P>();
            numAccept += 1;
        }
        else
            getAuxField()[freq, site] = save;
        numTotal += 1;
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R, P>(true);
        calcGreen<P>();

        numTotal = 0;
        numAccept = 0;
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_random(HamiltonMC<Tr>& hmc) {
        action.template randAuxField<R>();

        VectorND<Tr> init(getAuxField().getSize() * 2);
        init.read(reinterpret_cast<const Tr*>(getAuxField().data()));
        hmc.setInitial(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step(HamiltonMC<Tr>& hmc, bool warmup) {
        const auto& pos = hmc.template step<R, P>(*this);
        getAuxField().read(reinterpret_cast<const T*>(pos.data()));

        auto [lnAD, sgnD] = calcDet<P>();
        lnWeight = lnAD;
        sign = sgnD;
        if (!warmup)
            calcGreen<P>();
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step_for(int numStep, HamiltonMC<Tr>& hmc) {
        for (int i = 0; i < numStep; ++i)
            step<R, P>(hmc, true);
        calcGreen<P>();
        hmc.reset();
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::potentialV(const Vector auto& pos) -> Trv {
        assert(pos.getLength() == getAuxField().getSize() * 2 && "[Error]: Real matrix contains 2x number of elements of complex matrix");
        getAuxField().read(reinterpret_cast<const T*>(pos.data()));

        MatrixND<Tr> buffer = getAuxField().squaredNorms();
        buffer *= reciprocal(getBetaU());
        buffer.row(0) *= Trv(0.5);
        return buffer.sum() - calcDet<P>()[0];
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::potentialV(const Cell& cell) -> Trv {
        return potentialV<P>(cell.getPos().flatten());
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void FreqDQMC<T>::forceAsync(const Cell& cell, Vector auto& result) {
        assert(result.getLength() == getAuxField().getSize() * 2);
        getAuxField().read(reinterpret_cast<const T*>(cell.getPos().data()));
        for (auto& spinLU : lu) {
            spinLU.setWorking(action);
            action.flip();
        }

        Array<MatrixND<T>, 2> spinFs{};
        auto task = parallel_for<P>([this, &spinFs](size_t spin) {
            auto& spinLU = lu[spin];
            spinLU.compute();
            const MatrixND<T> inv = spinLU.inv();

            auto& spinF = spinFs[spin];
            spinF.resize(getAuxField());
            spinF.zeros();

            auto freqKernel = MatrixND<T>(2 * getNumFreq());
            Trv factor = spin == 0 ? 1.0 : -1.0;
            for (int site = 0; site < getNumSite(); ++site) {
                makeFreqKernel(freqKernel, inv, site);
                spinF[0, site].real() += freqKernel.diag().reals().sum() * factor;
                for (int delta = 1; delta < 2 * getNumFreq(); ++delta) {
                    T f1 = freqKernel.diag(delta).sum();
                    T f2 = freqKernel.diag(-delta).sum();
                    spinF[delta, site].real() += (f1.real() + f2.real()) * factor;
                    spinF[delta, site].imag() += (f1.imag() - f2.imag()) * factor;
                }
            }
        }, 2);

        MatrixND<T> force = -getAuxField() * (Trv(2) / getBetaU());
        force.row(0) *= Trv(0.5);
        task.wait();

        for (auto& spinF : spinFs)
            force += spinF;
        result.read(reinterpret_cast<const Tr*>(force.data()));
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::calcDet() -> Vector2D<Trv> {
        for (auto& spinLU : lu) {
            spinLU.setWorking(action);
            action.flip();
        }

        parallel_for<P>([this](size_t spin) {
            lu[spin].compute();
        }, 2).wait();

        Trv lnAD = 0, sgnD = 1;
        for (auto& spinLU : lu) {
            lnAD += ln(abs(spinLU.getMatrixLU().diag())).sum();
            sgnD *= unit(spinLU.getMatrixLU().diag()).prod().real();
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::calcGreen() {
        return parallel_for<P>([this](size_t spin) {
            auto& spinLU = lu[spin];
            MatrixND<T> inv = spinLU.inv();

            const int numSite = getNumSite();
            auto& green = greens[spin];
            green.zeros();
            for (int _ = 0, offset = 0; _ < action.getNumFreq() * 2; ++_) {
                green += inv.block(offset, numSite, offset, numSite).reals();
                offset += numSite;
            }
            green.diag() += Trv(0.5);
        }, 2);
    }

    template<Scalar T>
    int FreqDQMC<T>::calcFreqCutoff(Trv beta, Trv freqDensity) noexcept {
        int i = static_cast<int>((beta * freqDensity).toMachine() * 0.25);
        return std::max(i, 1) * 4; // Multiple of 4 so that SIMD works
    }

    template<Scalar T>
    void FreqDQMC<T>::makeFreqKernel(MatrixND<T>& freqKernel, const MatrixND<T>& invAction, int site) const {
        const int numSite = getNumSite();
        const int size = 2 * getNumFreq();
        assert(freqKernel.isSquare() && freqKernel.getRow() == size);
        for (int r = 0; r < size; ++r)
            for (int c = 0; c < size; ++c)
                freqKernel[r, c] = invAction.transpose().calc(r * numSite + site, c * numSite + site);
    }

    template<Scalar T>
    auto FreqDQMC<T>::getBetaU() const noexcept -> Trv {
        return getParams().getBeta() * getParams().getRepelU();;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::calcLnWeight() -> Vector2D<Trv> {
        Trv betaU = getBetaU();
        auto [lnAD, sgnD] = calcDet<P>();
        lnAD = lnAD
             - ln1pexp(lncosh(getAuxField().row(0).reals()) + fma(betaU, Trv(-0.5), MathConst<Trv>::ln2)).sum()
             - ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().reals()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum()
             - ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().imags()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum();
        return {lnAD, sgnD};
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<FreqDQMC<T>> {
    public:
        constexpr static bool IsPeriodBoundary = false;
    };
}
