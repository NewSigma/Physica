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

#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "DQMCImpl/ActionMatrix.h"
#include "ElasticDQMC.h"

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
        Trv correction;

        GreenPair greens;
        Trv lnWeight = Trv::nan();
        Trv sign = 1;
    public:
        FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity, int maxBoson);
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
        bool step(bool warmup = false);
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep);

        template<RNG R>
        void step_random(HamiltonMC<Tr>& hmc);
        template<RNG R, ExecutePolicy P = Sequential>
        Trv step(HamiltonMC<Tr>& hmc);
        [[nodiscard]] VectorND<Trv> makeDefaultMass() const;

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Vector auto& pos, Vector auto& result);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return action.getParams(); }
        [[nodiscard]] auto&& getAuxField(this auto&& self) noexcept { return self.action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return getParams().getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return action.getMaxBoson(); }
        [[nodiscard]] const auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Trv getSign() const noexcept { return sign; }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return false; }
    private:
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
    FreqDQMC<T>::FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity, int maxBoson)
            : action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity), maxBoson)
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
    }

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity))
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
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
    bool FreqDQMC<T>::step(bool warmup) {
        assert(lnWeight.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const int freq = std::uniform_int_distribution<int>(0, getMaxBoson() - 1)(R::getInstance());
        const T save = action.template randAuxField<R>(freq, site);

        auto [lnW, sgnD] = calcLnWeight<P>();
        const bool accept = Trv::template random_uniform<R>() < exp(lnW - lnWeight);
        if (accept) {
            lnWeight = lnW;
            sign = sgnD;
            if (!warmup)
                calcGreen<P>();
        }
        else
            getAuxField()[freq, site] = save;
        return accept;
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R, P>(true);
        calcGreen<P>();
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_random(HamiltonMC<Tr>& hmc) {
        action.template randAuxField<R>();

        VectorND<Tr> init(getAuxField().getSize() * 2);
        init.read(getAuxField());
        hmc.setInitPosition(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto FreqDQMC<T>::step(HamiltonMC<Tr>& hmc) -> Trv {
        Trv acceptR = hmc.template step<R, P>(*this);
        getAuxField().read(hmc.getSample());

        auto [lnAD, sgnD] = calcDet<P>();
        lnWeight = lnAD;
        sign = sgnD;
        calcGreen<P>();
        return acceptR;
    }

    template<Scalar T>
    auto FreqDQMC<T>::makeDefaultMass() const -> VectorND<Trv> {
        constexpr int Major = std::remove_cvref_t<decltype(getAuxField())>::getMajor(); // FIXME: clang 22 rejects getAuxField().getMajor()?
        VectorND<Trv> result(getAuxField().getSize() * 2, 1);
        auto mat = result.template reshape<Major>(getAuxField().getRow(), getAuxField().getCol() * 2);
        mat.bottomRows(1) *= Trv(2);
        return result;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto FreqDQMC<T>::potentialV(const Vector auto& pos) -> Trv {
        assert(pos.getLength() == getAuxField().getSize() * 2 && "[Error]: Real matrix contains 2x number of elements of complex matrix");
        getAuxField().read(pos);
        bool noAuxField = getBetaU().isSubNormal();
        if (noAuxField)
            return -calcDet<P>()[0];

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
    void FreqDQMC<T>::forceAsync(const Vector auto& pos, Vector auto& result) {
        assert(result.getLength() == pos.getLength());
        assert(result.getLength() == getAuxField().getSize() * 2);
        getAuxField().read(pos);
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

            const Trv factor = spin == 0 ? 1.0 : -1.0;
            const int size = getMaxBoson();
            const int numSite = getNumSite();
            for (int r = 0; r < size; ++r) {
                for (int c = 0; c < size; ++c) {
                    const auto diag = inv.transpose().block(r * numSite, numSite, c * numSite, numSite).diag();
                    if (r == c)
                        spinF.row(0) += diag.reals() * factor;
                    else if (r < c)
                        spinF.row(c - r) += diag * factor;
                    else
                        spinF.row(r - c) += diag.conjugate() * factor;
                }
            }
        }, 2);

        MatrixND<T> force(getAuxField().getRow(), getAuxField().getCol());
        bool noAuxField = getBetaU().isSubNormal();
        if (noAuxField)
            force.zeros();
        else {
            force = -getAuxField() * (Trv(2) / getBetaU());
            force.row(0) *= Trv(0.5);
        }
        task.wait();

        for (auto& spinF : spinFs)
            force += spinF;
        result.read(force);
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void FreqDQMC<T>::forceAsync(const Cell& cell, Vector auto& result) {
        forceAsync<P>(cell.getPos().flatten(), result);
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
    auto FreqDQMC<T>::calcLnWeight() -> Vector2D<Trv> {
        Trv betaU = getBetaU();
        auto [lnAD, sgnD] = calcDet<P>();
        lnAD = lnAD - ln1pexp(lncosh(getAuxField().row(0).reals()) + fma(betaU, Trv(-0.5), MathConst<Trv>::ln2)).sum();
        if (getMaxBoson() > 1) {
            lnAD -= ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().reals()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum()
                  + ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().imags()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum();
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
            for (int _ = 0, offset = 0; _ < 2 * getNumFreq(); ++_) {
                green += inv.block(offset, numSite, offset, numSite).reals();
                offset += numSite;
            }
            green.diag() += Trv(0.5) + correction;
        }, 2);
    }

    template<Scalar T>
    auto FreqDQMC<T>::getBetaU() const noexcept -> Trv {
        return getParams().getBeta() * getParams().getRepelU();
    }
}
