/*
 * Copyright 2023 WeiBo He.
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
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<class ScalarType, class PosScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class FreeModel {
        using PlainScalar = typename ScalarType::PlainScalar;
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim, NumReplica>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using MassVector = typename MDCellType::MassVector;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
        using BufferType = typename RingPolymerType::BufferType;
        using VectorType = Vector<ScalarType>;
        using Vector2D = Vector<ScalarType, 2>;
    private:
        VectorType omegaK;
        ScalarType omegaW;
        Utils::Array<Vector2D> coeffMatrixBase;
        ScalarType lastTimeStep;
    public:
        FreeModel();
        FreeModel(ScalarType temperatureT, size_t numReplica);
        FreeModel(const FreeModel&) = default;
        FreeModel(FreeModel&&) noexcept = default;
        ~FreeModel() = default;
        /* Operators */
        FreeModel& operator=(FreeModel obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        template<class Barostat>
        void npt_step(RingPolymerType& ringPolymer, MDCellType& cell, Barostat& barostat, ScalarType deltaT);
        void swap(FreeModel& obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getOmegaW() const noexcept { return omegaW; }
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

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::FreeModel() : lastTimeStep(0) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::FreeModel(ScalarType temperatureT, size_t numReplica)
            : omegaW(RingPolymerType::calcOmegaW(temperatureT, numReplica))
            , lastTimeStep(0) {
        const size_t kSpaceSize = RingPolymerType::calcKSpaceSize(numReplica);
        omegaK.resize(kSpaceSize);
        coeffMatrixBase.resize(kSpaceSize);
        for (size_t i = 0; i < omegaK.getLength(); ++i)
            omegaK[i] = omegaW * sin(PlainScalar(M_PI * i / numReplica)) * PlainScalar(2);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
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

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    template<class Barostat>
    void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::npt_step(
            RingPolymerType& ringPolymer,
            MDCellType& cell,
            Barostat& barostat,
            ScalarType deltaT) {
        nve_step(ringPolymer, deltaT);
        LatticeMatrix lattice = cell->getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        cell.setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    FreeModel<ScalarType, PosScalarType, Dim, NumReplica>&
    FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::operator=(FreeModel<ScalarType, PosScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::swap(FreeModel& obj) noexcept {
        omegaK.swap(obj.omegaK);
        omegaW.swap(obj.omegaW);
        coeffMatrixBase.swap(obj.coeffMatrixBase);
        lastTimeStep.swap(obj.lastTimeStep);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    inline void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::pre_nve_step_impl(
            [[maybe_unused]] RingPolymerType& ringPolymer, ScalarType deltaT) {
        assert(NumReplica != 1);
        assert(omegaK.getLength() == ringPolymer.getBuffer().getColumn());
        if (lastTimeStep != deltaT)
            updateTimeStep(deltaT);
    }
    /**
     * Integrate strategy as introduced in [1]
     * 
     * Reference:
     * [1] M. Ceriotti, M. Parrinello, T. E. Markland and D. E. Manolopoulos, J. Chem. Phys. 133, 124104 (2010).
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    inline void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::do_nve_step_impl(
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
            const ScalarType factor = mass * omegaK[j];
            const ScalarType cosine = coeffMatrixBase[j][0];
            const ScalarType sine = coeffMatrixBase[j][1];
            auto col = buffer.col(j);
            const auto momentum = col[0];
            const auto pos = col[1];
            const auto new_momentum = cosine * momentum - (factor * sine) * pos;
            const auto new_pos = (sine / factor) * momentum + cosine * pos;
            col[0] = new_momentum;
            col[1] = new_pos;
        }
        ringPolymer.toBeadRepr(id_dof, output);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::nve_step_impl(
            RingPolymerType& ringPolymer, const PhaseMatrix& input, PhaseMatrix& output, ScalarType deltaT) {
        pre_nve_step_impl(ringPolymer, deltaT);

        const size_t dof = input.getRow() / 2U;
        for (size_t i = 0; i < dof; ++i)
            do_nve_step_impl(i, ringPolymer, input, output);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica>
    void FreeModel<ScalarType, PosScalarType, Dim, NumReplica>::updateTimeStep(ScalarType deltaT) {
        for (size_t i = 0; i < omegaK.getLength(); ++i) {
            const ScalarType phase = omegaK[i] * deltaT;
            coeffMatrixBase[i] = Vector2D{cos(phase), sin(phase)};
        }
        lastTimeStep = deltaT;
    }
}
