/*
 * Copyright 2023-2024 Weibo He.
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
     * Exact as introduced in [1]: Exact free ring-polymer update
     * Cayley as introduced in [2]: More accurate and efficient Cayley transform-based approximated version of BAOAB
     * 
     * Reference:
     * [1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
     * [2] J. Chem. Phys. 152, 104102 (2020); https://doi.org/10.1063/1.5134810
     */
    enum class RPMDIntegrator {
        Exact,
        Cayley
    };

    template<class ScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<class ScalarType,
             unsigned int Dim,
             size_t NumReplica,
             RPMDIntegrator Integrator>
    class FreeModel {
        using PlainScalar = typename ScalarType::PlainScalar;
        using ComplexType = ComplexScalar<ScalarType>;
        using MDType = RPMD<ScalarType, Dim, NumReplica>;
        using MDCellType = typename MDType::MDCellType;
        using RingPolymerType = typename MDType::RingPolymerType;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using MassVector = typename MDCellType::MassVector;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
        using BufferType = typename RingPolymerType::BufferType;
        using FFTType = typename RingPolymerType::FFTType;
        using VectorType = Vector<ScalarType>;
        using Vector2D = Vector<ScalarType, 2>;
        static_assert(std::is_same<ComplexType, typename BufferType::ScalarType>::value, "[Error]: Inconsistent type");
    private:
        VectorType omegaK;
        ScalarType omegaW;
        Array<Vector2D> coeffMatrixBase;
        ScalarType lastTimeStep;
        size_t numReplica;
    public:
        FreeModel();
        FreeModel(ScalarType temperatureT, size_t numReplica_);
        FreeModel(const FreeModel&) = default;
        FreeModel(FreeModel&&) noexcept = default;
        ~FreeModel() = default;
        /* Operators */
        FreeModel& operator=(FreeModel obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        template<class Barostat>
        void npt_step(MDType& rpmd, Barostat& barostat, ScalarType deltaT);
        void swap(FreeModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getOmegaW() const noexcept { return omegaW; }
        /* Setters */
        void setTemperature(ScalarType temperatureT);
    protected:
        inline void pre_nve_step_impl([[maybe_unused]] RingPolymerType& ringPolymer, ScalarType deltaT);
        inline void do_nve_step_impl(
                size_t id_dof,
                RingPolymerType& ringPolymer,
                const PhaseMatrix& input,
                PhaseMatrix& output);
        void nve_step_impl(RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output, ScalarType deltaT);
    private:
        void updateTimeStep(ScalarType deltaT);
    };

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    FreeModel<ScalarType, Dim, NumReplica, Integrator>::FreeModel() : lastTimeStep(0) {}

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    FreeModel<ScalarType, Dim, NumReplica, Integrator>::FreeModel(ScalarType temperatureT, size_t numReplica_)
            : lastTimeStep(0)
            , numReplica(numReplica_) {
        const size_t kSpaceSize = FFTType::rSizeToKSize(numReplica);
        omegaK.resize(kSpaceSize);
        coeffMatrixBase.resize(kSpaceSize);
        setTemperature(temperatureT);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
        if constexpr (NumReplica != 1)
            nve_step_impl(ringPolymer, ringPolymer.asMatrix(), ringPolymer.asMatrix(), deltaT);
        else {
            const size_t dof = ringPolymer.getDOF();
            const auto& massVec = ringPolymer.getMassVec();
            auto& phase = ringPolymer.asMatrix();
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                phase(i + dof, 0) += phase(i, 0) * deltaT / mass;
            }
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    template<class Barostat>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::npt_step(
            MDType& rpmd,
            Barostat& barostat,
            ScalarType deltaT) {
        nve_step(rpmd.getRingPolymer(), deltaT);
        LatticeMatrix lattice = rpmd.getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        rpmd.setLattice(std::move(lattice));
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    FreeModel<ScalarType, Dim, NumReplica, Integrator>&
    FreeModel<ScalarType, Dim, NumReplica, Integrator>::operator=(
            FreeModel<ScalarType, Dim, NumReplica, Integrator> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::swap(FreeModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        omegaK.swap(obj.omegaK);
        omegaW.swap(obj.omegaW);
        coeffMatrixBase.swap(obj.coeffMatrixBase);
        lastTimeStep.swap(obj.lastTimeStep);
        std::swap(numReplica, obj.numReplica);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::setTemperature(ScalarType temperatureT) {
        omegaW = RingPolymerType::calcOmegaW(temperatureT, numReplica);
        for (size_t i = 0; i < omegaK.getLength(); ++i)
            omegaK[i] = omegaW * sin(PlainScalar(M_PI * i / numReplica)) * PlainScalar(2);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    inline void FreeModel<ScalarType, Dim, NumReplica, Integrator>::pre_nve_step_impl(
            [[maybe_unused]] RingPolymerType& ringPolymer, ScalarType deltaT) {
        assert(NumReplica != 1);
        assert(omegaK.getLength() == ringPolymer.getBuffer().getColumn());
        if (lastTimeStep != deltaT)
            updateTimeStep(deltaT);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    inline void FreeModel<ScalarType, Dim, NumReplica, Integrator>::do_nve_step_impl(
            size_t id_dof,
            RingPolymerType& ringPolymer,
            const PhaseMatrix& input,
            PhaseMatrix& output) {
        const auto& massVec = ringPolymer.getMassVec();
        const size_t kSpaceSize = omegaK.getLength();
        BufferType& buffer = ringPolymer.getBuffer();

        const auto mass = massVec[id_dof / Dim];
        ringPolymer.toNormalRepr(id_dof, input);
        /* Translational mode */ {
            buffer(1, 0) = buffer(1, 0) + buffer(0, 0) * (lastTimeStep / mass);
        }
        for (size_t j = 1; j < kSpaceSize; ++j) {
            auto col = buffer.col(j);
            const ComplexType momentum = col[0];
            const ComplexType pos = col[1];
            ComplexType new_momentum, new_pos;
            if constexpr (Integrator == RPMDIntegrator::Exact) {
                const ScalarType factor = mass * omegaK[j];
                const ScalarType cosine = coeffMatrixBase[j][0];
                const ScalarType sine = coeffMatrixBase[j][1];
                new_momentum = cosine * momentum - (factor * sine) * pos;
                new_pos = (sine / factor) * momentum + cosine * pos;
            }
            else {
                const ScalarType factor = coeffMatrixBase[j][0];
                const ScalarType factorDeltaT = coeffMatrixBase[j][1];
                new_momentum = factor * momentum - (mass * factorDeltaT * square(omegaK[j])) * pos;
                new_pos = (factorDeltaT / mass) * momentum + factor * pos;
            }
            col[0] = new_momentum;
            col[1] = new_pos;
        }
        ringPolymer.toBeadRepr(id_dof, output);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::nve_step_impl(
            RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output, ScalarType deltaT) {
        pre_nve_step_impl(ringPolymer, deltaT);

        const size_t dof = input.getRow() / 2U;
        for (size_t i = 0; i < dof; ++i)
            do_nve_step_impl(i, ringPolymer, input, output);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, RPMDIntegrator Integrator>
    void FreeModel<ScalarType, Dim, NumReplica, Integrator>::updateTimeStep(ScalarType deltaT) {
        for (size_t i = 0; i < omegaK.getLength(); ++i) {
            const ScalarType phase = omegaK[i] * deltaT;
            if constexpr (Integrator == RPMDIntegrator::Exact)
                coeffMatrixBase[i] = Vector2D{cos(phase), sin(phase)};
            else {
                const ScalarType factor = reciprocal(sqrt(ScalarType(1) + square(phase)));
                coeffMatrixBase[i] = Vector2D{factor, factor * deltaT};
            }
        }
        lastTimeStep = deltaT;
    }
}

namespace Physica {
    template<class T1, unsigned int T2, size_t T3, RPMDIntegrator Integrator>
    class Traits<Core::FreeModel<T1, T2, T3, Integrator>> {
    public:
        using ScalarType = T1;
        constexpr static unsigned int Dim = T2;
        constexpr static size_t NumReplica = T3;
        constexpr static bool IsPeriodBoundary = true;
    };
}
