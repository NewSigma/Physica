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
    template<class ScalarType, class PosScalarType, unsigned int Dim> class RPMD;

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class DynamicStep {
        using RPMDType = RPMD<ScalarType, PosScalarType, Dim>;
        using LatticeMatrix = typename RPMDType::LatticeMatrix;
        using BufferType = typename RPMDType::BufferType;
        using VectorType = Vector<ScalarType>;
        using Vector2D = Vector<ScalarType, 2>;

        RPMDType* rpmd;
        VectorType omegaK;
        Utils::Array<Vector2D> coeffMatrixBase;
        ScalarType lastTimeStep;
    public:
        DynamicStep();
        DynamicStep(RPMDType& rpmd_);
        DynamicStep(const DynamicStep&) = default;
        DynamicStep(DynamicStep&&) noexcept = default;
        ~DynamicStep() = default;
        /* Operators */
        DynamicStep& operator=(DynamicStep obj) noexcept;
        /* Operations */
        void nve_step(ScalarType deltaT);
        template<class Barostat>
        void npt_step(Barostat& barostat, ScalarType deltaT);
        void swap(DynamicStep& obj) noexcept;
    private:
        void updateTimeStep(ScalarType deltaT);
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>::DynamicStep() : rpmd(nullptr), lastTimeStep(0) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>::DynamicStep(RPMDType& rpmd_)
            : rpmd(&rpmd_)
            , omegaK(rpmd_.buffer.getColumn())
            , coeffMatrixBase(rpmd_.buffer.getColumn())
            , lastTimeStep(0) {
        const size_t numReplica = rpmd->getNumReplica();
        for (size_t i = 0; i < omegaK.getLength(); ++i)
            omegaK[i] = rpmd->getOmegaW() * sin(ScalarType(M_PI * i / numReplica)) * 2;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void DynamicStep<ScalarType, PosScalarType, Dim>::nve_step(ScalarType deltaT) {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using ComplexVector2D = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = rpmd->getDOF();
        const size_t kSpaceSize = omegaK.getLength();
        if (lastTimeStep != deltaT)
            updateTimeStep(deltaT);

        MatrixType matA(2, 2);
        ComplexVector2D temp{};
        BufferType& buffer = rpmd->buffer;
        assert(kSpaceSize == buffer.getColumn());
        
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = rpmd->getMassVec()[i / Dim];
            rpmd->toNormalRepr(i);
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
            rpmd->toBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Barostat>
    void DynamicStep<ScalarType, PosScalarType, Dim>::npt_step(
            Barostat& barostat,
            ScalarType deltaT) {
        nve_step(rpmd, deltaT);
        LatticeMatrix lattice = rpmd->getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        rpmd->setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    DynamicStep<ScalarType, PosScalarType, Dim>&
    DynamicStep<ScalarType, PosScalarType, Dim>::operator=(DynamicStep<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void DynamicStep<ScalarType, PosScalarType, Dim>::swap(DynamicStep& obj) noexcept {
        std::swap(rpmd, obj.rpmd);
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
