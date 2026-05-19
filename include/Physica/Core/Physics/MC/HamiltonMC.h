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

#include "Physica/Core/Physics/MD/KineticModel/OpenModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica {
    /**
     * Implements the NUTS sampler introduced in [1]
     *
     * Reference:
     * [1] J. Mach. Learn. Res. 15(1), 1593–1623 (2014); https://dl.acm.org/doi/10.5555/2627435.2638586
     */
    template<Scalar T>
    class HamiltonMC {
        using This = HamiltonMC<T>;
        using Node = RPMD<T, 1, 1>;
        using KineticModel = OpenModel<T, 1, 1>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;
        struct Proposal {
            VectorND<T> curX;
            int numAccept;
            bool converge;
            Trv acceptR; // Averaged acceptance rate during traverse
            int numVisited;
        };

        Node root;
        Node nodeF;
        Node nodeR;
        KineticModel kinetic;
        Trv targetR;
        Trv maxDelta;

        VectorND<T> sample;
        Trv upperE;
    public:
        HamiltonMC(VectorND<T> mass, Trv targetAcceptRate = 0.65, Trv maxDelta = 1000); // Default value from [1]
        HamiltonMC(const This&) = default;
        HamiltonMC(This&&) noexcept = default;
        ~HamiltonMC() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<RNG R, ExecutePolicy P = Sequential>
        VectorND<Trv> warmup(int numWarmup, auto&& forceModel);
        template<RNG R, ExecutePolicy P = Sequential>
        Trv step(auto&& forceModel);
        /* Getters */
        [[nodiscard]] auto& getRoot() noexcept { return root; }
        [[nodiscard]] size_t getDOF() const noexcept { return root.getDOF(); }
        [[nodiscard]] Trv getTargetR() const noexcept { return targetR; }
        [[nodiscard]] const auto& getSample() const noexcept { return sample; }
        /* Setters */
        void setInitPosition(VectorND<T> init) noexcept;
        void setTimeStep(Trv timestep) noexcept;
    private:
        template<RNG R, ExecutePolicy P = Sequential>
        void initTimeStep(auto&& forceModel);
        template<RNG R, ExecutePolicy P>
        Proposal traverse(bool forward, int height, auto&& forceModel);

        bool canProgress(const Node& forward, const Node& reverse) const noexcept;
        /* Static members */
        static MDCell<T, 1> makeDummyCell(VectorND<T> mass);
    };

    template<Scalar T>
    HamiltonMC<T>::HamiltonMC(VectorND<T> mass, Trv targetAcceptRate, Trv maxDelta)
            : root(makeDummyCell(std::move(mass)), 1, 1, 1, 1)
            , nodeF(root)
            , nodeR(root)
            , kinetic(1, 1)
            , targetR(targetAcceptRate)
            , maxDelta(maxDelta)
            , sample(getDOF(), 0) {
        assert(maxDelta.isPositive());
        assert(targetAcceptRate.isPositive() && targetAcceptRate <= Trv(1));
        root.getPhaseMatrix().zeros();
    }
    /**
     * Algo. 6 of [1]
     */
    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto HamiltonMC<T>::warmup(int numWarmup, auto&& forceModel) -> VectorND<Trv> {
        initTimeStep<R, P>(forceModel);

        VectorND<Trv> lnTimeSteps(numWarmup);
        Trv initialMu = ln(Trv(10) * root.getTimeStep());
        Trv lnConjTimeStep = 0;
        Trv weightAcceptR = 0;
        for (int i = 1; i <= numWarmup; ++i) {
            constexpr Trv Gamma = 0.05;
            constexpr Trv Kappa = 0.75;
            constexpr int Shift = 10;
            Trv acceptR = step<R, P>(forceModel);
            Trv factor = reciprocal(Trv(i + Shift));
            weightAcceptR = fma(factor, targetR - acceptR, (Trv(1) - factor) * weightAcceptR);

            Trv m = Trv(i);
            Trv lnTimeStep = initialMu - weightAcceptR * sqrt(m) / Gamma;
            root.setTimeStep(exp(lnTimeStep));

            factor = pow(m, -Kappa);
            lnConjTimeStep = fma(factor, lnTimeStep, (Trv(1) - factor) * lnConjTimeStep);
            lnTimeSteps[i - 1] = lnConjTimeStep;
        }
        root.setTimeStep(exp(lnConjTimeStep));
        return lnTimeSteps;
    }
    /**
     * Algo. 3 of [1]
     */
    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto HamiltonMC<T>::step(auto&& forceModel) -> Trv {
        sample.assign(root.getPhaseMatrix().col(0).tail(getDOF()));
        root.template initMomentum<KineticModel, R>();
        root.template updateForce<P>(forceModel);
        nodeF = root;
        nodeR = root;
        upperE = root.template calcClassicalInternalEnergy<P>(forceModel)
               - ln(Trv::template random_uniform<R>() + std::numeric_limits<T>::min());

        int iteration = 0;
        int numAccept = 1;
        while (true) {
            bool forward = R::coin();
            auto [curX, numAccept_, converge, acceptR, _] = traverse<R, P>(forward, iteration, forceModel);
            if (!converge && (numAccept_ > 0)) {
                Trv prob = Trv(numAccept_) / Trv(numAccept);
                bool accept = Trv::template random_uniform<R>() < prob;
                if (accept)
                    sample = std::move(curX);
            }

            if (converge || !canProgress(nodeF, nodeR))
                return acceptR;

            numAccept += numAccept_;
            iteration += 1;
        }
    }

    template<Scalar T>
    void HamiltonMC<T>::setInitPosition(VectorND<T> init) noexcept {
        assert(getDOF() == init.getLength());
        sample = std::move(init);
    }

    template<Scalar T>
    void HamiltonMC<T>::setTimeStep(Trv timestep) noexcept {
        root.setTimeStep(timestep);
    }
    /**
     * Algo. 4 of [1]
     */
    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void HamiltonMC<T>::initTimeStep(auto&& forceModel) {
        root.template initMomentum<KineticModel, R>();
        const Trv e0 = root.template calcClassicalInternalEnergy<P>(forceModel);

        Trv prev = 0;
        Trv timestep = 1;
        while (true) {
            root.setTimeStep(timestep);
            nodeF = root;
            nodeF.template nve_step<P>(kinetic, forceModel);
            Trv e = nodeF.template calcClassicalInternalEnergy<P>(forceModel);
            Trv factor = (exp(e0 - e) > 0.5) ? Trv(2) : Trv(0.5);
            bool converge = prev.isPositive() && (prev != factor);
            if (converge)
                break;

            prev = factor;
            timestep *= factor;
        }
    }
    /**
     * Traverse the balanced tree starting from \param node at \param height; moving in the \param forward direction.
     *
     * \returns a proposed sample with metadata; one may accept or reject it
     */
    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto HamiltonMC<T>::traverse(bool forward, int height, auto&& forceModel) -> Proposal {
        auto& node = forward ? nodeF : nodeR;
        if (height == 0) { // Leaf node
            Trv prevE = node.template calcClassicalInternalEnergy<P>(forceModel);
            if (forward)
                node.template nve_step<P>(kinetic, forceModel);
            else
                node.template nve_step_back<P>(kinetic, forceModel);
            root = node;

            Trv curE = node.template calcClassicalInternalEnergy<P>(forceModel);
            return {
                .curX = node.getPhaseMatrix().col(0).tail(getDOF()),
                .numAccept = curE < upperE,
                .converge = curE >= upperE + maxDelta,
                .acceptR = std::min(exp(prevE - curE), Trv(1)),
                .numVisited = 1
            };
        }

        assert(height > 0);
        auto [curX, numAccept, converge, acceptR, numVisited] = traverse<R, P>(forward, height - 1, forceModel);
        if (!converge) {
            auto [curX_, numAccept_, converge_, acceptR_, numVisited_] = traverse<R, P>(forward, height - 1, forceModel);
            if (numAccept_ > 0) {
                Trv prob = Trv(numAccept_) / Trv(numAccept + numAccept_);
                bool accept = Trv::template random_uniform<R>() < prob;
                if (accept)
                    curX = std::move(curX_);
            }

            bool process = forward ? canProgress(node, root) : canProgress(root, node);
            numAccept += numAccept_;
            converge = converge_ || !process;
            acceptR = [=]() -> Trv {
                Trv total = numVisited + numVisited_;
                Trv factor1 = Trv(numVisited) / total;
                Trv factor2 = Trv(numVisited_) / total;
                return fma(acceptR, factor1, acceptR_ * factor2);
            }();
            numVisited += numVisited_;
        }
        return {std::move(curX), numAccept, converge, acceptR, numVisited};
    }

    template<Scalar T>
    MDCell<T, 1> HamiltonMC<T>::makeDummyCell(VectorND<T> mass) {
        size_t dof = mass.getLength();
        return MDCell<T, 1>(dof, std::move(mass));
    }

    template<Scalar T>
    bool HamiltonMC<T>::canProgress(const Node& forward, const Node& reverse) const noexcept {
        auto phaseF = forward.getPhaseMatrix().col(0);
        auto phaseR = reverse.getPhaseMatrix().col(0);
        auto momentF = phaseF.head(getDOF());
        auto momentR = phaseR.head(getDOF());
        VectorND<T> delta = phaseF.tail(getDOF()) - phaseR.tail(getDOF());
        return (delta * momentF).isPositive() && (delta * momentR).isPositive();
    }
}
