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
    template<class ScalarType, class PosScalarType, unsigned int Dim> class RPMD;

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class DynamicStep {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using MassVector = typename MDCellType::MassVector;
        using BufferType = typename RingPolymerType::BufferType;
        using VectorType = Vector<ScalarType>;
        using Vector2D = Vector<ScalarType, 2>;

        VectorType omegaK;
        Utils::Array<Vector2D> coeffMatrixBase;
        ScalarType lastTimeStep;
    public:
        DynamicStep();
        DynamicStep(ScalarType omegaW, size_t kSpaceSize, size_t numReplica);
        DynamicStep(const DynamicStep&) = default;
        DynamicStep(DynamicStep&&) noexcept = default;
        ~DynamicStep() = default;
        /* Operators */
        DynamicStep& operator=(DynamicStep obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, const MassVector& mass, ScalarType deltaT);
        template<class Barostat>
        void npt_step(RingPolymerType& ringPolymer, MDCellType& cell, Barostat& barostat, ScalarType deltaT);
        void swap(DynamicStep& obj) noexcept;
    private:
        void updateTimeStep(ScalarType deltaT);
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>::DynamicStep() : lastTimeStep(0) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>::DynamicStep(ScalarType omegaW, size_t kSpaceSize, size_t numReplica)
            : omegaK(kSpaceSize)
            , coeffMatrixBase(kSpaceSize)
            , lastTimeStep(0) {
        for (size_t i = 0; i < omegaK.getLength(); ++i)
            omegaK[i] = omegaW * sin(ScalarType(M_PI * i / numReplica)) * 2;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void DynamicStep<ScalarType, PosScalarType, Dim>::nve_step(RingPolymerType& ringPolymer, const MassVector& massVec, ScalarType deltaT) {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using ComplexVector2D = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = ringPolymer.getDOF();
        const size_t kSpaceSize = omegaK.getLength();
        if (lastTimeStep != deltaT)
            updateTimeStep(deltaT);

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

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Barostat>
    void DynamicStep<ScalarType, PosScalarType, Dim>::npt_step(
            RingPolymerType& ringPolymer,
            MDCellType& cell,
            Barostat& barostat,
            ScalarType deltaT) {
        nve_step(ringPolymer, cell.getMassVec(), deltaT);
        LatticeMatrix lattice = cell->getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        cell.setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>&
    DynamicStep<ScalarType, PosScalarType, Dim>::operator=(DynamicStep<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void DynamicStep<ScalarType, PosScalarType, Dim>::swap(DynamicStep& obj) noexcept {
        omegaK.swap(obj.omegaK);
        coeffMatrixBase.swap(obj.coeffMatrixBase);
        lastTimeStep.swap(obj.lastTimeStep);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void DynamicStep<ScalarType, PosScalarType, Dim>::updateTimeStep(ScalarType deltaT) {
        for (size_t i = 0; i < omegaK.getLength(); ++i) {
            const ScalarType phase = omegaK[i] * deltaT;
            coeffMatrixBase[i] = Vector2D{cos(phase), sin(phase)};
        }
        lastTimeStep = deltaT;
    }
}
