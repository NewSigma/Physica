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
            VectorND<T> sample;
            Trv acceptR; // Averaged acceptance rate during visit
            int numAccept;
            int numVisited;
            bool stop; // If trajectory is not good
        };

        Node root;
        Node nodeF;
        Node nodeR;
        std::stack<Node> nodes;
        KineticModel kinetic;
        Trv targetR;
        Trv maxDelta;
        int maxTreeDepth;

        VectorND<T> sample;
        Trv upperE;
    public:
        HamiltonMC(VectorND<T> mass, Trv targetAcceptRate = 0.65, Trv maxDelta = 1000, int maxTreeDepth = 10); // Default value from [1]
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
        void initTimeStep(auto& forceModel);
        template<RNG R, ExecutePolicy P>
        [[nodiscard]] Proposal visit(int height, auto& forceModel);

        template<ExecutePolicy P>
        [[nodiscard]] Proposal visitLeaf(bool forward, auto& forceModel);
        /* Static members */
        [[nodiscard]] static MDCell<T, 1> makeDummyCell(VectorND<T> mass);
        template<RNG R>
        static void metropolis(Proposal& lhs, const Proposal& rhs, bool process) noexcept;
        static bool canProgress(const Node& forward, const Node& reverse) noexcept;
    };

    template<Scalar T>
    HamiltonMC<T>::HamiltonMC(VectorND<T> mass, Trv targetAcceptRate, Trv maxDelta, int maxTreeDepth)
            : root(makeDummyCell(std::move(mass)), 1, 1, 1, 1)
            , nodeF(root)
            , nodeR(root)
            , kinetic(1, 1)
            , targetR(targetAcceptRate)
            , maxDelta(maxDelta)
            , maxTreeDepth(maxTreeDepth)
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
            auto [sample_, acceptR, numAccept_, _, stop] = visit<R, P>(iteration, forceModel);
            if (!stop && (numAccept_ > 0)) {
                Trv prob = Trv(numAccept_) / Trv(numAccept);
                bool accept = Trv::template random_uniform<R>() < prob;
                if (accept)
                    sample = std::move(sample_);
            }

            if (stop || !canProgress(nodeF, nodeR))
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
    void HamiltonMC<T>::initTimeStep(auto& forceModel) {
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
     * Visit the balanced tree of height \param height; moving in the \param forward direction.
     *
     * \returns a proposed sample with metadata; one may accept or reject it
     */
    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto HamiltonMC<T>::visit(int height, auto& forceModel) -> Proposal {
        return [this, forward = R::coin(), height, &forceModel](this const auto& self) noexcept -> Proposal {
            bool isLeaf = height == nodes.size();
            if (isLeaf)
                return visitLeaf<P>(forward, forceModel);

            assert(height > 0);
            nodes.emplace(forward ? nodeF : nodeR);
            auto lhs = self();
            if (!lhs.stop) {
                const auto rhs = self();
                const auto& saved = nodes.top();
                bool process = forward ? canProgress(nodeF, saved) : canProgress(saved, nodeR);
                metropolis<R>(lhs, rhs, process);
            }
            nodes.pop();
            return lhs;
        }();
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto HamiltonMC<T>::visitLeaf(bool forward, auto& forceModel) -> Proposal {
        auto& node = forward ? nodeF : nodeR;
        Trv prevE = node.template calcClassicalInternalEnergy<P>(forceModel);
        if (forward)
            node.template nve_step<P>(kinetic, forceModel);
        else
            node.template nve_step_back<P>(kinetic, forceModel);

        Trv curE = node.template calcClassicalInternalEnergy<P>(forceModel);
        return Proposal{
            .sample = node.getPhaseMatrix().col(0).tail(getDOF()),
            .acceptR = std::min(exp(prevE - curE), Trv(1)),
            .numAccept = curE < upperE,
            .numVisited = 1,
            .stop = (curE >= upperE + maxDelta) || (nodes.size() >= maxTreeDepth)};
    }

    template<Scalar T>
    MDCell<T, 1> HamiltonMC<T>::makeDummyCell(VectorND<T> mass) {
        size_t dof = mass.getLength();
        return MDCell<T, 1>(dof, std::move(mass));
    }

    template<Scalar T>
    template<RNG R>
    void HamiltonMC<T>::metropolis(Proposal& lhs, const Proposal& rhs, bool process) noexcept {
        auto& [sample1, acceptR1, numAccept1, numVisited1, stop1] = lhs;
        const auto& [sample2, acceptR2, numAccept2, numVisited2, stop2] = rhs;
        if (numAccept2 > 0) {
            Trv prob = Trv(numAccept2) / Trv(numAccept1 + numAccept2);
            bool accept = Trv::template random_uniform<R>() < prob;
            if (accept)
                sample1 = std::move(sample2);
        }
        numAccept1 += numAccept2;
        stop1 = stop2 || !process;
        acceptR1 = fma(acceptR1, Trv(numVisited1), acceptR2 * Trv(numVisited2)) / Trv(numVisited1 + numVisited2);
        numVisited1 += numVisited2;
    }

    template<Scalar T>
    bool HamiltonMC<T>::canProgress(const Node& forward, const Node& reverse) noexcept {
        assert(forward.getDOF() == reverse.getDOF());
        size_t dof = forward.getDOF();
        auto phaseF = forward.getPhaseMatrix().col(0);
        auto phaseR = reverse.getPhaseMatrix().col(0);
        auto momentF = phaseF.head(dof);
        auto momentR = phaseR.head(dof);
        VectorND<T> delta = phaseF.tail(dof) - phaseR.tail(dof);
        return (delta * momentF).isPositive() && (delta * momentR).isPositive();
    }
}
