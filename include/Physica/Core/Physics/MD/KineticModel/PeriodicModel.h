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
    template<class ScalarType, class PosScalarType, unsigned int Dim> class RingPolymer;

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class PeriodicModel {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using MassVector = typename MDCellType::MassVector;
        using BufferType = typename RingPolymerType::BufferType;
        using VectorType = Vector<ScalarType>;
        using Vector2D = Vector<ScalarType, 2>;

        VectorType omegaK;
        ScalarType omegaW;
        Utils::Array<Vector2D> coeffMatrixBase;
        ScalarType lastTimeStep;
    public:
        PeriodicModel();
        PeriodicModel(ScalarType temperatureT, size_t numReplica);
        PeriodicModel(const PeriodicModel&) = default;
        PeriodicModel(PeriodicModel&&) noexcept = default;
        ~PeriodicModel() = default;
        /* Operators */
        PeriodicModel& operator=(PeriodicModel obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        template<class Barostat>
        void npt_step(RingPolymerType& ringPolymer, MDCellType& cell, Barostat& barostat, ScalarType deltaT);
        void swap(PeriodicModel& obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getOmegaW() const noexcept { return omegaW; }
    private:
        void updateTimeStep(ScalarType deltaT);
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    PeriodicModel<ScalarType, PosScalarType, Dim>::PeriodicModel() : lastTimeStep(0) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    PeriodicModel<ScalarType, PosScalarType, Dim>::PeriodicModel(ScalarType temperatureT, size_t numReplica)
            : omegaW(RingPolymerType::calcOmegaW(temperatureT, numReplica))
            , lastTimeStep(0) {
        const size_t kSpaceSize = RingPolymerType::calcKSpaceSize(numReplica);
        omegaK.resize(kSpaceSize);
        coeffMatrixBase.resize(kSpaceSize);
        for (size_t i = 0; i < omegaK.getLength(); ++i)
            omegaK[i] = omegaW * sin(ScalarType(M_PI * i / numReplica)) * 2;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void PeriodicModel<ScalarType, PosScalarType, Dim>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using ComplexVector2D = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = ringPolymer.getDOF();
        const size_t kSpaceSize = omegaK.getLength();
        const auto& massVec = ringPolymer.getMassVec();
        if (lastTimeStep != deltaT)
            updateTimeStep(deltaT);

        if (ringPolymer.getNumReplica() > 1) {
            MatrixType matA(2, 2);
            ComplexVector2D temp{};
            BufferType& buffer = ringPolymer.getBuffer();
            assert(kSpaceSize == buffer.getColumn());
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                ringPolymer.toNormalRepr(i);
                /* Translational mode */ {
                    buffer(1, 0) += buffer(0, 0) * deltaT / mass;
                }
                for (size_t j = 1; j < kSpaceSize; ++j) {
                    auto col = buffer.col(j);
                    const ScalarType factor = ScalarType(mass) * omegaK[j];
                    const ScalarType cosine = coeffMatrixBase[j][0];
                    const ScalarType sine = coeffMatrixBase[j][1];
                    matA(0, 0) = cosine;
                    matA(0, 1) = -factor * sine;
                    matA(1, 0) = sine / factor;
                    matA(1, 1) = cosine;
                    temp = matA * col;
                    col = temp;
                }
                ringPolymer.toBeadRepr(i);
            }
        }
        else {
            auto& phase = ringPolymer.asMatrix();
            for (size_t i = 0; i < dof; ++i) {
                const auto mass = massVec[i / Dim];
                phase(i + dof, 0) += phase(i, 0) * deltaT / mass;
            }
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Barostat>
    void PeriodicModel<ScalarType, PosScalarType, Dim>::npt_step(
            RingPolymerType& ringPolymer,
            MDCellType& cell,
            Barostat& barostat,
            ScalarType deltaT) {
        nve_step(ringPolymer, deltaT);
        LatticeMatrix lattice = cell->getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        cell.setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    PeriodicModel<ScalarType, PosScalarType, Dim>&
    PeriodicModel<ScalarType, PosScalarType, Dim>::operator=(PeriodicModel<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void PeriodicModel<ScalarType, PosScalarType, Dim>::swap(PeriodicModel& obj) noexcept {
        omegaK.swap(obj.omegaK);
        omegaW.swap(obj.omegaW);
        coeffMatrixBase.swap(obj.coeffMatrixBase);
        lastTimeStep.swap(obj.lastTimeStep);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void PeriodicModel<ScalarType, PosScalarType, Dim>::updateTimeStep(ScalarType deltaT) {
        for (size_t i = 0; i < omegaK.getLength(); ++i) {
            const ScalarType phase = omegaK[i] * deltaT;
            coeffMatrixBase[i] = Vector2D{cos(phase), sin(phase)};
        }
        lastTimeStep = deltaT;
    }
}
