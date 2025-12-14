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

#include "DQMCImpl/ActionMatrix.h"
#include "DQMCImpl/ImagKinetic.h"

namespace Physica {
    template<Scalar T>
    class FreqDQMC {
        using This = FreqDQMC<T>;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        using GreenPair = ImagKinetic<Tr>::GreenPair;
    private:
        Array<DenseLU<T, false>, 2> lu;
        ActionMatrix<T> action;
        MatrixND<T> solBuffer;

        GreenPair greens;
        Trv lnAbsDet = Trv::nan();
        Trv sign = 1;
    public:
        FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity);
        FreqDQMC(const This&) = default;
        FreqDQMC(This&&) noexcept = default;
        ~FreqDQMC() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R>
        void step_random();

        template<RNG R>
        void step();
        template<RNG R>
        void step_for(int numStep);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return action.getParams(); }
        [[nodiscard]] int getNumSite() const noexcept { return getParams().getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Trv getSign() const noexcept { return sign; }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        /* Static members */
        [[nodiscard]] static int calcFreqCutoff(Trv beta, Trv freqDensity) noexcept;
    private:
        [[nodiscard]] Vector2D<Trv> calcDet();
        void calcGreen();
    };

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, calcFreqCutoff(params.getBeta(), freqDensity))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        solBuffer.resize(size);
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_random() {
        action.template randAuxField<R>();
        auto [lnAD, sgnD] = calcDet();
        lnAbsDet = lnAD;
        sign = sgnD;
        calcGreen();
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step() {
        assert(lnAbsDet.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const int freq = std::uniform_int_distribution<int>(0, getNumFreq() * 2 - 1)(R::getInstance());
        const T save = action.template randAuxField<R>(freq, site);

        auto [lnAD, sgnD] = calcDet();
        const bool accept = Trv::template random_uniform<R>() < exp(lnAD - lnAbsDet);
        if (accept) {
            lnAbsDet = lnAD;
            sign = sgnD;
            calcGreen();
        }
        else
            action.getAuxField()(freq, site) = save;
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R>();
    }

    template<Scalar T>
    int FreqDQMC<T>::calcFreqCutoff(Trv beta, Trv freqDensity) noexcept {
        int i = static_cast<int>((beta * freqDensity).toMachine() * 0.25);
        return std::max(i, 1) * 4; // Multiple of 4 so that SIMD works
    }

    template<Scalar T>
    auto FreqDQMC<T>::calcDet() -> Vector2D<Trv> {
        lu[0].compute(action);

        action.flip();
        lu[1].compute(action);

        Trv lnAD = 0, sgnD = 1;
        for (const auto& spinLU : lu) {
            lnAD += ln(abs(spinLU.getMatrixLU().diag())).sum();
            sgnD *= unit(spinLU.getMatrixLU().diag()).prod().real();
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    void FreqDQMC<T>::calcGreen() {
        const int numSite = getNumSite();
        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            solBuffer = spinLU.inv();

            auto& green = greens[spin];
            green.zeros();
            for (int _ = 0, offset = 0; _ < action.getNumFreq() * 2; ++_) {
                green += solBuffer.block(offset, numSite, offset, numSite).reals();
                offset += numSite;
            }
            green *= reciprocal(getParams().getBeta());
            green.diag() += Trv(0.5);
        }
    }
}