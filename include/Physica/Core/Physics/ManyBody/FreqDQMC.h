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
        HubbardParams<Tr> params;
        Array<DenseLU<T, false>, 2> lu;
        ActionMatrix<T> action;
        HamiltonMC<Tr> hmc;
        Trv correction;

        GreenPair greens;
        Trv lnWeight = Trv::nan();
        Trv sign = 1;
    public:
        FreqDQMC(HubbardParams<Tr> params_, Trv freqDensity, int maxBoson);
        FreqDQMC(HubbardParams<Tr> params_, Trv freqDensity);
        FreqDQMC(const This&) = default;
        FreqDQMC(This&&) noexcept = default;
        ~FreqDQMC() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R>
        void step_random();
        template<RNG R, ExecutePolicy P = Sequential>
        Trv step();
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep);

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Vector auto& pos, Vector auto& result);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);

        Tr calcBerry();
        /* Getters */
        [[nodiscard]] auto&& getAuxField(this auto&& self) noexcept { return self.action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return action.getMaxBoson(); }
        [[nodiscard]] auto&& getHMC(this auto&& self) noexcept { return self.hmc; }
        [[nodiscard]] const auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Trv getSign() const noexcept { return sign; }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return false; }
    private:
        template<ExecutePolicy P>
        [[nodiscard]] Vector2D<Trv> calcDet();
        template<ExecutePolicy P>
        auto calcGreen();
        /* Getters */
        [[nodiscard]] Trv getBetaU() const noexcept;
        /* Static members */
        [[nodiscard]] static VectorND<Trv> makeDefaultMass(const Matrix auto& auxField);
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(HubbardParams<Tr> params_, Trv freqDensity, int maxBoson)
            : params(std::move(params_))
            , action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity), maxBoson)
            , hmc(makeDefaultMass(getAuxField()))
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
    }

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(HubbardParams<Tr> params_, Trv freqDensity)
            : params(std::move(params_))
            , action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity))
            , hmc(makeDefaultMass(getAuxField()))
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_random() {
        action.template random_normal<R>();

        VectorND<Tr> init(getAuxField().getSize() * 2);
        init.read(getAuxField());
        hmc.setInitPosition(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto FreqDQMC<T>::step() -> Trv {
        Trv acceptR = hmc.template step<R, P>(*this);
        bool hasAuxField = !getBetaU().isSubNormal();
        if (hasAuxField)
            getAuxField().read(hmc.getSample());
        else
            getAuxField().zeros();

        auto [lnAD, sgnD] = calcDet<P>();
        lnWeight = lnAD;
        sign = sgnD;
        calcGreen<P>();
        return acceptR;
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void FreqDQMC<T>::step_for(int numStep) {
        for (int _ = 0; _ < numStep; ++_)
            step<R, P>();
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

            const Trv factor = spin == 0 ? 1 : -1;
            const int numFreq2 = getNumFreq() * 2;
            for (int site = 0; site < getNumSite(); ++site) {
                const auto block = inv.transpose().block(site * numFreq2, numFreq2, site * numFreq2, numFreq2);
                for (int freq = 0; freq < getMaxBoson(); ++freq) {
                    if (freq == 0)
                        spinF[freq, site] = block.diag().sum().real() * factor;
                    else
                        spinF[freq, site] = (block.diag(freq).sum() + block.diag(-freq).conjugate().sum()) * factor;
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
    auto FreqDQMC<T>::calcBerry() -> Tr {
        const int numSite = getNumSite();
        const int numFreq2 = 2 * getNumFreq();
        MatrixND<T> inv(numSite * numFreq2);
        VectorND<Tr> ptrace(numFreq2);
        Tr result = 0;
        for (auto& spinLU : lu) {
            inv = spinLU.inv();

            ptrace.zeros();
            for (int site = 0; site < numSite; ++site)
                ptrace += inv.block(site * numFreq2, numFreq2, site * numFreq2, numFreq2).diag().imags();
            result -= ptrace * action.getMatsubara().diag();
        }
        const Tr background = reciprocal(Trv(1) + square(divide(params.calcBetaMu(), action.getMatsubara().diag().head(getNumFreq())))).sum() * (4 * numSite);
        return result - background;
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
            lnAD += spinLU.lnAbsDet();
            sgnD *= spinLU.sgndet().real();
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
        return params.getBeta() * params.getRepelU();
    }

    template<Scalar T>
    auto FreqDQMC<T>::makeDefaultMass(const Matrix auto& auxField) -> VectorND<Trv> {
        constexpr int Major = auxField.getMajor();
        VectorND<Trv> result(auxField.getSize() * 2, 1);
        auto mat = result.template reshape<Major>(auxField.getRow(), auxField.getCol() * 2);
        mat.bottomRows(1) *= Trv(2);
        return result;
    }
}
