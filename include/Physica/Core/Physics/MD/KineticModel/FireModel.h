/*
 * Copyright 2023 Weibo He.
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

namespace Physica::Core {
    /**
     * Fast inertial relaxation engine(FIRE) as introduced in [1, 2]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 97, 170201 (2006); https://doi.org/10.1103/PhysRevLett.97.170201
     * [2] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584
     */
    template<Scalar T, unsigned int Dim>
    class FireModel {
    public:
        using MDType = RPMD<T, Dim, 1>;
    private:
        T timeStep;
        T maxTimeStep;
        T stepIncFactor;
        T stepDecFactor;
        T mixAlpha0;
        T mixAlpha;
        T alphaDecFactor;
        unsigned int numLasyStep;

        T normF;
        unsigned int numStep;
    public:
        FireModel(T timeStep_,
                  T maxTimeStep,
                  T stepIncFactor_ = 1.1,
                  T stepDecFactor_ = 0.5,
                  T mixAlpha_ = 0.1,
                  T alphaDecFactor_ = 0.99,
                  unsigned int numLasyStep_ = 5);
        FireModel(const FireModel&) = default;
        FireModel(FireModel&&) noexcept = default;
        ~FireModel() = default;
        /* Operators */
        FireModel& operator=(FireModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void paramStep(MDType& rpmd);
        void mixingStep(MDType& rpmd);
        void swap(FireModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getTimeStep() const noexcept { return timeStep; }
        [[nodiscard]] T getMixAlpha() const noexcept { return mixAlpha; }
        [[nodiscard]] T getForceNorm() const noexcept { return normF; } //Use force norm as convergence criteria as suggested by [2]
    };
    /**
     * Default param provided in [1]
     */
    template<Scalar T, unsigned int Dim>
    FireModel<T, Dim>::FireModel(
            T timeStep_,
            T maxTimeStep_,
            T stepIncFactor_,
            T stepDecFactor_,
            T mixAlpha_,
            T alphaDecFactor_,
            unsigned int numLasyStep_)
            : timeStep(timeStep_)
            , maxTimeStep(maxTimeStep_)
            , stepIncFactor(stepIncFactor_)
            , stepDecFactor(stepDecFactor_)
            , mixAlpha0(mixAlpha_)
            , alphaDecFactor(alphaDecFactor_)
            , numLasyStep(numLasyStep_)
            , numStep(0) {
        assert(timeStep.isPositive());
        assert(timeStep <= maxTimeStep);
        assert(stepIncFactor > T(1));
        assert(T(0) < stepDecFactor && stepDecFactor < T(1));
        assert(numLasyStep > 0 && "[Error]: It is unstable if numLasyStep = 0");
        assert(T(0) < mixAlpha0 && mixAlpha0 < T(1));
        assert(T(0) < alphaDecFactor && alphaDecFactor < T(1));
        mixAlpha = mixAlpha0;
    }

    template<Scalar T, unsigned int Dim>
    void FireModel<T, Dim>::paramStep(MDType& rpmd) {
        const auto force = rpmd.getForce().col(0);
        auto phase = rpmd.getPhaseMatrix().col(0);
        auto momentum = phase.head(rpmd.getDOF());
        const T power = force * momentum; //Assuming unit mass
        if (power.isPositive()) {
            if (numStep > numLasyStep) {
                timeStep = std::min(stepIncFactor * timeStep, maxTimeStep);
                mixAlpha *= alphaDecFactor;
            }
            numStep += 1;
        }
        else {
            timeStep *= stepDecFactor;
            momentum = T(0);
            mixAlpha = mixAlpha0;
            numStep = 0;
        }
    }

    template<Scalar T, unsigned int Dim>
    void FireModel<T, Dim>::mixingStep(MDType& rpmd) {
        const auto force = rpmd.getForce().col(0);
        auto phase = rpmd.getPhaseMatrix().col(0);
        auto momentum = phase.head(rpmd.getDOF());

        const T normP = momentum.norm();
        normF = force.norm();
        momentum = (T(1) - mixAlpha) * momentum + (mixAlpha * normP / normF) * force;
    }

    template<Scalar T, unsigned int Dim>
    void FireModel<T, Dim>::swap(FireModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        timeStep.swap(obj.timeStep);
        maxTimeStep.swap(obj.maxTimeStep);
        stepIncFactor.swap(obj.stepIncFactor);
        stepDecFactor.swap(obj.stepDecFactor);
        mixAlpha0.swap(obj.mixAlpha0);
        mixAlpha.swap(obj.mixAlpha);
        alphaDecFactor.swap(obj.alphaDecFactor);
        std::swap(numLasyStep, obj.numLasyStep);

        normF.swap(obj.normF);
        std::swap(numStep, obj.numStep);
    }
}
