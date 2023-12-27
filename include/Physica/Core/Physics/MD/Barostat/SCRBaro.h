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

#include "Berendsen.h"
#include "Physica/Core/Math/Calculus/ODE/SRK2.h"

namespace Physica::Core {
    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type> class SCRBaro;

    namespace Internal {
        template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
        class Traits<SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>> {
        public:
            constexpr static unsigned int Order = 1;
        };
    }
    /**
     * Stochastic cell rescaling(SCR) barostat as introduced in [1], its anisotropic version as introduced in [2]
     * 
     * Reference:
     * [1] J. Chem. Phys. 153, 114107 (2020); https://doi.org/10.1063/5.0020514
     * [2] https://doi.org/10.48550/arXiv.2111.06402
     */
    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    class SCRBaro : private Berendsen<ScalarType, NumReplica, Type> {
        using Base = Berendsen<ScalarType, NumReplica, Type>;
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::Vector3D;
        using Base::Dim;
    public:
        SCRBaro() = default;
        SCRBaro(ScalarType compressRate, ScalarType tempT, ScalarType targetP);
        SCRBaro(const SCRBaro&) = default;
        SCRBaro(SCRBaro&&) noexcept = default;
        ~SCRBaro() = default;
        /* Operators */
        SCRBaro& operator=(SCRBaro obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class MDType, class ForceModel>
        void npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT);
        using Base::swap;
        /* Getters */
        using Base::getLastStress;
        /* Setters */
        using Base::setTemperature;
    private:
        [[nodiscard]] LatticeMatrix makeDiffuseMatrix(ScalarType pressPerDOF) const;
        [[nodiscard]] LatticeMatrix makeDeltaLattice(const LatticeMatrix& lattice, ScalarType volume, ScalarType deltaT) const;
    };

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::SCRBaro(
            ScalarType compressRate, ScalarType tempT, ScalarType targetP) : Base(compressRate, tempT, targetP) {}

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    template<class MDType, class ForceModel>
    void SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::npt_step(
            MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT) {
        Base::lastStress = stress;

        const size_t numReplica = rpmd.getNumReplica();
        const size_t numParticle = rpmd.getNumParticle();
        const size_t dof = rpmd.getDOF();
        auto& phase = rpmd.getPhaseMatrix();
        const auto deltaLattice = makeDeltaLattice(rpmd.getLattice(), rpmd.getVolume(), deltaT);
        const LatticeMatrix scaleMatrix = Base::makeScaleMatrix(rpmd.getInvLattice(), deltaLattice);
        for (size_t i = 0; i < numReplica; ++i) {
            auto col = phase.col(i);
            auto momentum = col.head(dof);
            for (size_t j = 0; j < numParticle; ++j) {
                auto momentum_j = momentum.segment(j * Dim, (j + 1) * Dim);
                const Vector3D delta = scaleMatrix * momentum_j;
                momentum_j -= delta;
            }
            auto pos = col.tail(dof);
            for (size_t j = 0; j < numParticle; ++j) {
                auto pos_j = pos.segment(j * Dim, (j + 1) * Dim);
                const Vector3D delta = scaleMatrix * pos_j;
                pos_j += delta;
            }
        }
        rpmd.setLattice(rpmd.getLattice() + deltaLattice);
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDiffuseMatrix(ScalarType pressPerDOF) const {
        const ScalarType diffuseFactor = sqrt(ScalarType(2.0 / Dim) * Base::compressRate * pressPerDOF);
        auto& gen = RandomPoolType::getGen();
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (Type == BaroType::Anisotropic) {
            const auto rand = LatticeMatrix::random_normal(Dim, Dim, gen);
            result = (rand + rand.transpose()) * (diffuseFactor * ScalarType(0.5));
        }
        else if constexpr (Type == BaroType::XY) {
            result(0, 0) = ScalarType::random_normal(gen);
            result(1, 1) = ScalarType::random_normal(gen);
            result(0, 1) = result(1, 0) = ScalarType::random_normal(gen) * ScalarType(M_SQRT1_2);
            auto corner = result.topLeftCorner(2);
            corner *= diffuseFactor;
        }
        else if constexpr (Type == BaroType::Z) {
            result(2, 2) = ScalarType::random_normal(gen) * diffuseFactor;
        }
        else {
            constexpr bool Unreachable = Type == BaroType::Anisotropic;
            static_assert(Unreachable, "[Error]: Not implemented");
        }
        return result;
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDeltaLattice(
            const LatticeMatrix& lattice, ScalarType volume, ScalarType deltaT) const {
        const ScalarType pressPerDOF = (Base::tempT * PhyConst<AU>::boltzmannK) / volume;
        const auto decayMatrix = Base::makeDecayMatrix(pressPerDOF);
        const auto diffuseMatrix = makeDiffuseMatrix(pressPerDOF);
        LatticeMatrix result(Dim, Dim, 0);
        auto integrateKernel = [&, deltaT](size_t r, size_t c) {
            using Integrator = SRK2<ScalarType, 1>;
            using VectorType = typename Integrator::VectorType;
            [[maybe_unused]] ScalarType unused = 0;
            VectorType sol{ScalarType(0)};
            Integrator::step(deltaT, unused, sol,
                [&lattice, &decayMatrix, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] VectorType sol) -> VectorType {
                    return {(decayMatrix * lattice).calc(r, c)};
                },
                [&lattice, &diffuseMatrix, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] VectorType sol) -> VectorType {
                    return {(diffuseMatrix * lattice).calc(r, c)};
                });
            result(r, c) = sol[0];
        };

        if constexpr (Type == BaroType::Anisotropic) {
            for (size_t r = 0; r < Dim; ++r)
                for (size_t c = 0; c < Dim; ++c)
                    integrateKernel(r, c);
        }
        else if constexpr (Type == BaroType::XY) {
            for (size_t r = 0; r < 2; ++r)
                for (size_t c = 0; c < 2; ++c)
                    integrateKernel(r, c);
        }
        else if constexpr (Type == BaroType::Z) {
            integrateKernel(2, 2);
        }
        else {
            constexpr bool Unreachable = Type == BaroType::Anisotropic;
            static_assert(Unreachable, "[Error]: Not implemented");
        }
        return result;
    }
}
