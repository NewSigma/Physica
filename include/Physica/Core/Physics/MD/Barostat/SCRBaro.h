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

#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Math/Calculus/ODE/SRK2.h"
#include "BaroType.h"

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
    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    class SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic> {
        constexpr static unsigned int Dim = 3;
        using MDType = RPMD<ScalarType, Dim, NumReplica>;
        using LatticeMatrix = typename MDType::MDCellType::LatticeMatrix;

        ScalarType compressRate;
        ScalarType tempT;
        ScalarType targetP;
        LatticeMatrix lastStress;
    public:
        SCRBaro() = default;
        SCRBaro(ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_);
        SCRBaro(const SCRBaro&) = default;
        SCRBaro(SCRBaro&&) noexcept = default;
        ~SCRBaro() = default;
        /* Operators */
        SCRBaro& operator=(SCRBaro obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class ForceModel>
        void npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT);
        void swap(SCRBaro& obj) noexcept;
    private:
        [[nodiscard]] LatticeMatrix makeDecayMatrix(ScalarType pressPerDOF) const;
        [[nodiscard]] LatticeMatrix makeDiffuseMatrix(ScalarType pressPerDOF) const;
    };

    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::SCRBaro(
            ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_)
            : compressRate(compressRate_)
            , tempT(tempT_)
            , targetP(targetP_) {
        assert(compressRate.isPositive() && "[Error]: Compress rate must be positive");
        assert(tempT.isPositive() && "[Error]: Temperature must be positive");
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    template<class ForceModel>
    void SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::npt_step(
            MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT) {
        lastStress = stress;
        const ScalarType pressPerDOF = (tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = makeDecayMatrix(pressPerDOF);
        const auto diffuseMatrix = makeDiffuseMatrix(pressPerDOF);
        const auto& lattice = rpmd.getLattice();
        LatticeMatrix deltaLattice(Dim, Dim);
        for (size_t r = 0; r < Dim; ++r) {
            for (size_t c = 0; c < Dim; ++c) {
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
                deltaLattice(r, c) = sol[0];
            }
        }
        //Scale pos and momentum
        const size_t numReplica = rpmd.getNumReplica();
        const size_t numParticle = rpmd.getNumParticle();
        const size_t dof = rpmd.getDOF();
        auto& phase = rpmd.getPhaseMatrix();
        const LatticeMatrix scaleMat = LatticeMatrix::unitMatrix(Dim) + deltaLattice * rpmd.getInvLattice();
        const LatticeMatrix invScaleMat = LatticeMatrix::unitMatrix(Dim) - deltaLattice * rpmd.getInvLattice();
        for (size_t i = 0; i < numReplica; ++i) {
            auto col = phase.col(i);
            auto momentum = col.head(dof);
            for (size_t j = 0; j < numParticle; ++j) {
                auto momentum_j = momentum.segment(j * Dim, (j + 1) * Dim);
                momentum_j = invScaleMat.transpose() * momentum_j;
            }
            auto pos = col.tail(dof);
            for (size_t j = 0; j < numParticle; ++j) {
                auto pos_j = pos.segment(j * Dim, (j + 1) * Dim);
                pos_j = scaleMat * pos_j;
            }
        }
        rpmd.setLattice(lattice + deltaLattice);
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    void SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::swap(SCRBaro& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        compressRate.swap(obj.compressRate);
        tempT.swap(obj.tempT);
        targetP.swap(obj.targetP);
        lastStress.swap(obj.lastStress);
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::makeDecayMatrix(ScalarType pressPerDOF) const {
        LatticeMatrix result = lastStress;
        const ScalarType centerTargetP = targetP - pressPerDOF;
        for (size_t i = 0; i < Dim; ++i)
            result(i, i) -= centerTargetP;
        result *= compressRate / ScalarType(Dim);
        return result;
    }
    
    template<class ScalarType, size_t NumReplica, class RandomPoolType>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, BaroType::Anisotropic>::makeDiffuseMatrix(ScalarType pressPerDOF) const {
        LatticeMatrix result;
        const auto rand = LatticeMatrix::random_normal(Dim, Dim, RandomPoolType::getGen());
        const ScalarType diffuseFactor = sqrt(ScalarType(2.0 / Dim) * compressRate * pressPerDOF);
        result = (rand + rand.transpose()) * (diffuseFactor * ScalarType(0.5));
        return result;
    }
}
