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
    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    class SCRBaro {
        constexpr static unsigned int Dim = 3;
        using MDType = RPMD<ScalarType, Dim, NumReplica>;
        using LatticeMatrix = typename MDType::MDCellType::LatticeMatrix;
        using Vector3D = Vector<ScalarType, Dim>;

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
        void swap(SCRBaro& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLastStress() const noexcept { return lastStress; }
    private:
        [[nodiscard]] LatticeMatrix makeDecayMatrix(ScalarType pressPerDOF) const;
        [[nodiscard]] LatticeMatrix makeDiffuseMatrix(ScalarType pressPerDOF) const;
        [[nodiscard]] LatticeMatrix makeDeltaLattice(MDType& rpmd, ScalarType deltaT) const;
        [[nodiscard]] LatticeMatrix makeScaleMatrix(const MDType& rpmd, const LatticeMatrix& deltaLattice) const;
    };

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::SCRBaro(
            ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_)
            : compressRate(compressRate_)
            , tempT(tempT_)
            , targetP(targetP_) {
        assert(compressRate.isPositive() && "[Error]: Compress rate must be positive");
        assert(tempT.isPositive() && "[Error]: Temperature must be positive");
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    template<class ForceModel>
    void SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::npt_step(
            MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT) {
        lastStress = stress;

        const size_t numReplica = rpmd.getNumReplica();
        const size_t numParticle = rpmd.getNumParticle();
        const size_t dof = rpmd.getDOF();
        auto& phase = rpmd.getPhaseMatrix();
        const auto deltaLattice = makeDeltaLattice(rpmd, deltaT);
        const LatticeMatrix scaleMatrix = makeScaleMatrix(rpmd, deltaLattice);
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
    void SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::swap(SCRBaro& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        compressRate.swap(obj.compressRate);
        tempT.swap(obj.tempT);
        targetP.swap(obj.targetP);
        lastStress.swap(obj.lastStress);
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDecayMatrix(ScalarType pressPerDOF) const {
        LatticeMatrix result = lastStress;
        const ScalarType centerTargetP = targetP - pressPerDOF;
        if constexpr (Type == BaroType::Anisotropic) {
            for (size_t i = 0; i < Dim; ++i)
                result(i, i) -= centerTargetP;
        }
        else {
            for (size_t i = 0; i < 2; ++i)
                result(i, i) -= centerTargetP;
            auto col = result.col(2);
            col = ScalarType(0);
            auto row = result.row(2);
            row = ScalarType(0);
        }
        result *= compressRate / ScalarType(Dim);
        return result;
    }
    
    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDiffuseMatrix(ScalarType pressPerDOF) const {
        const ScalarType diffuseFactor = sqrt(ScalarType(2.0 / Dim) * compressRate * pressPerDOF);
        auto& gen = RandomPoolType::getGen();
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (Type == BaroType::Anisotropic) {
            const auto rand = LatticeMatrix::random_normal(Dim, Dim, gen);
            result = (rand + rand.transpose()) * (diffuseFactor * ScalarType(0.5));
        }
        else {
            result(0, 0) = ScalarType::random_normal(gen);
            result(1, 1) = ScalarType::random_normal(gen);
            result(0, 1) = result(1, 0) = ScalarType::random_normal(gen) + ScalarType::random_normal(gen);
            auto corner = result.topLeftCorner(2);
            corner *= diffuseFactor;
        }
        return result;
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDeltaLattice(MDType& rpmd, ScalarType deltaT) const {
        const ScalarType pressPerDOF = (tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = makeDecayMatrix(pressPerDOF);
        const auto diffuseMatrix = makeDiffuseMatrix(pressPerDOF);
        const auto& lattice = rpmd.getLattice();
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
        else {
            for (size_t r = 0; r < 2; ++r)
                for (size_t c = 0; c < 2; ++c)
                    integrateKernel(r, c);
        }
        return result;
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeScaleMatrix(
            const MDType& rpmd, const LatticeMatrix& deltaLattice) const {
        LatticeMatrix result;
        if constexpr (Type == BaroType::Anisotropic)
            result = deltaLattice * rpmd.getInvLattice();
        else {
            result = deltaLattice * rpmd.getInvLattice();
            auto col = result.col(2);
            col = ScalarType(0);
        }
        return result;
    }
}
