/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica {
    /**
     * Reference:
     * [1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
     * [2] J. Chem. Phys. 152, 104102 (2020); https://doi.org/10.1063/1.5134810
     */
    enum class RPMDIntegrator : char {
        Exact, // [1]: Exact update for free ring-polymer system
        Cayley, // [2]: Accurate and efficient Cayley-transform-based approximated version of BAOAB
    };
    /**
     * TODO: Separate non-classical part out
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator = RPMDIntegrator::Exact>
    class FreeModel {
        using This = FreeModel<T, Dim, NumReplica, Integrator>;
        using MDType = RPMD<T, Dim, NumReplica>;
        using RingPolymerType = MDType::RingPolymerType;
        using PhaseMatrix = RingPolymerType::PhaseMatrix;

        using Tv = T::ValueType;
        using Tc = Complex<T>;
        constexpr static bool IsClassical = NumReplica == 1;
    private:
        Array<Vector2D<T>> coeffMatrixBase;
        VectorND<T> omegaK;
        T omegaW;
        T lastTimeStep = 0;
        size_t replica = 1;
    public:
        FreeModel();
        FreeModel(T temperatureT, size_t numReplica);
        FreeModel(const This&) = default;
        FreeModel(This&&) noexcept = default;
        ~FreeModel() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, T deltaT);
        void npt_step(MDType& rpmd, auto& barostat, T deltaT);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getOmegaW() const noexcept;
        [[nodiscard]] size_t getNumReplica() const noexcept;
        /* Setters */
        void setTemperature(T temperatureT);
    protected:
        void nve_step_impl(RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output, T deltaT);
        void pre_nve_step_impl([[maybe_unused]] RingPolymerType& ringPolymer, T deltaT);
        void do_nve_step_impl(size_t id_dof, RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output);
    private:
        void updateCoeff(T deltaT);
    };

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    FreeModel<T, Dim, NumReplica, Integrator>::FreeModel(T temperatureT, size_t numReplica)
            : replica(numReplica) {
        assert(NumReplica == numReplica || NumReplica == Dynamic);
        if constexpr (!IsClassical) {
            using FFTType = RingPolymerType::FFTType;
            const size_t kSpaceSize = FFTType::rSizeToKSize(getNumReplica());
            coeffMatrixBase.resize(kSpaceSize);
            omegaK.resize(kSpaceSize);
        }
        setTemperature(temperatureT);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::nve_step(RingPolymerType& ringPolymer, T deltaT) {
        if constexpr (NumReplica != 1)
            nve_step_impl(ringPolymer, ringPolymer.asMatrix(), ringPolymer.asMatrix(), deltaT);
        else {
            const size_t dof = ringPolymer.getDOF();
            const auto& massVec = ringPolymer.getMassVec();
            auto& phase = ringPolymer.asMatrix();
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                phase[i + dof, 0] = fma(phase[i, 0] / mass, deltaT, phase[i + dof, 0]);
            }
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::npt_step(MDType& rpmd, auto& barostat, T deltaT) {
        using LatticeMatrix = MDType::LatticeMatrix;
        nve_step(rpmd.getRingPolymer(), deltaT);
        LatticeMatrix lattice = rpmd.getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        rpmd.setLattice(std::move(lattice));
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        omegaK.swap(obj.omegaK);
        omegaW.swap(obj.omegaW);
        coeffMatrixBase.swap(obj.coeffMatrixBase);
        lastTimeStep.swap(obj.lastTimeStep);
        std::swap(replica, obj.replica);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    T FreeModel<T, Dim, NumReplica, Integrator>::getOmegaW() const noexcept {
        static_assert(!IsClassical);
        return omegaW;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    size_t FreeModel<T, Dim, NumReplica, Integrator>::getNumReplica() const noexcept {
        if constexpr (NumReplica != Dynamic)
            return NumReplica;
        else
            return replica;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::setTemperature(T temperatureT) {
        if constexpr (!IsClassical) {
            omegaW = RingPolymerType::calcOmegaW(temperatureT, getNumReplica());
            for (size_t i = 0; i < omegaK.getLength(); ++i)
                omegaK[i] = omegaW * sin(Tv(MathConst<T>::pi * i / getNumReplica())) * Tv(2);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::nve_step_impl(
            RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output, T deltaT) {
        pre_nve_step_impl(ringPolymer, deltaT);

        const size_t dof = input.getRow() / 2U;
        for (size_t i = 0; i < dof; ++i)
            do_nve_step_impl(i, ringPolymer, input, output);
    }


    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::pre_nve_step_impl(
            [[maybe_unused]] RingPolymerType& ringPolymer, T deltaT) {
        static_assert(!IsClassical);
        assert(omegaK.getLength() == ringPolymer.getBuffer().getCol());
        if (lastTimeStep != deltaT) [[unlikely]]
            updateCoeff(deltaT);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::do_nve_step_impl(size_t id_dof, RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output) {
        using BufferType = RingPolymerType::BufferType;
        static_assert(std::is_same<Tc, typename BufferType::ScalarType>::value, "[Error]: Inconsistent type");
        const auto mass = ringPolymer.getMassVec()[id_dof / Dim];
        BufferType& buffer = ringPolymer.getBuffer();
        ringPolymer.toNormalRepr(id_dof, input);
        /* Translational mode */ {
            buffer[1, 0] = buffer[1, 0] + buffer[0, 0] * (lastTimeStep / mass);
        }

        const size_t kSpaceSize = omegaK.getLength();
        for (size_t j = 1; j < kSpaceSize; ++j) {
            auto col = buffer.col(j);
            const Tc momentum = col[0];
            const Tc pos = col[1];
            Tc new_momentum, new_pos;
            if constexpr (Integrator == RPMDIntegrator::Exact) {
                const T factor = mass * omegaK[j];
                const T cosine = coeffMatrixBase[j][0];
                const T sine = coeffMatrixBase[j][1];
                new_momentum = cosine * momentum - (factor * sine) * pos;
                new_pos = (sine / factor) * momentum + cosine * pos;
            }
            else {
                const T factor = coeffMatrixBase[j][0];
                const T factorDeltaT = coeffMatrixBase[j][1];
                new_momentum = factor * momentum - (mass * factorDeltaT * square(omegaK[j])) * pos;
                new_pos = (factorDeltaT / mass) * momentum + factor * pos;
            }
            col[0] = new_momentum;
            col[1] = new_pos;
        }
        ringPolymer.toBeadRepr(id_dof, output);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<T, Dim, NumReplica, Integrator>::updateCoeff(T deltaT) {
        for (size_t i = 0; i < omegaK.getLength(); ++i) {
            const T phase = omegaK[i] * deltaT;
            if constexpr (Integrator == RPMDIntegrator::Exact)
                coeffMatrixBase[i] = Vector2D<T>{cos(phase), sin(phase)};
            else {
                const T factor = reciprocal(sqrt(T(1) + square(phase)));
                coeffMatrixBase[i] = Vector2D<T>{factor, factor * deltaT};
            }
        }
        lastTimeStep = deltaT;
    }
}

namespace Physica {
    template<Scalar T1, unsigned int T2, size_t T3, RPMDIntegrator Integrator>
    class Traits<FreeModel<T1, T2, T3, Integrator>> {
    public:
        using ScalarType = T1;
        constexpr static unsigned int Dim = T2;
        constexpr static size_t NumReplica = T3;
        constexpr static bool IsPeriodBoundary = true;
    };
}
