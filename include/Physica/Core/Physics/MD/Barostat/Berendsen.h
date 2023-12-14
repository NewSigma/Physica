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
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "BaroType.h"

namespace Physica::Core {
    template<class ScalarType, size_t NumReplica, BaroType Type> class Berendsen;

    namespace Internal {
        template<class ScalarType, size_t NumReplica, BaroType Type>
        class Traits<Berendsen<ScalarType, NumReplica, Type>> {
        public:
            constexpr static unsigned int Order = 1;
        };
    }
    /**
     * Berendsen barostat as introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 81, 3684 (1984); https://doi.org/10.1063/1.448118
     */
    template<class ScalarType, size_t NumReplica, BaroType Type>
    class Berendsen {
    public:
        constexpr static unsigned int Dim = 3;
        using MDType = RPMD<ScalarType, Dim, NumReplica>;
        using LatticeMatrix = typename MDType::MDCellType::LatticeMatrix;
        using Vector3D = Vector<ScalarType, Dim>;
    protected:
        ScalarType compressRate;
        ScalarType tempT;
        ScalarType targetP;
        LatticeMatrix lastStress;
    public:
        Berendsen() = default;
        Berendsen(ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_);
        Berendsen(const Berendsen&) = default;
        Berendsen(Berendsen&&) noexcept = default;
        ~Berendsen() = default;
        /* Operators */
        Berendsen& operator=(Berendsen obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class ForceModel>
        void npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT);
        void swap(Berendsen& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLastStress() const noexcept { return lastStress; }
        /* Setters */
        void setTemperature(ScalarType tempT_) { tempT = std::move(tempT_); }
    protected:
        [[nodiscard]] LatticeMatrix makeDecayMatrix(ScalarType pressPerDOF) const;
        [[nodiscard]] LatticeMatrix makeScaleMatrix(const MDType& rpmd, const LatticeMatrix& deltaLattice) const;
    private:
        [[nodiscard]] LatticeMatrix makeDeltaLattice(MDType& rpmd, ScalarType deltaT) const;
    };

    template<class ScalarType, size_t NumReplica, BaroType Type>
    Berendsen<ScalarType, NumReplica, Type>::Berendsen(
            ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_)
            : compressRate(compressRate_)
            , tempT(tempT_)
            , targetP(targetP_) {
        assert(compressRate.isPositive() && "[Error]: Compress rate must be positive");
        assert(tempT.isPositive() && "[Error]: Temperature must be positive");
    }

    template<class ScalarType, size_t NumReplica, BaroType Type>
    template<class ForceModel>
    void Berendsen<ScalarType, NumReplica, Type>::npt_step(
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

    template<class ScalarType, size_t NumReplica, BaroType Type>
    void Berendsen<ScalarType, NumReplica, Type>::swap(Berendsen& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        compressRate.swap(obj.compressRate);
        tempT.swap(obj.tempT);
        targetP.swap(obj.targetP);
        lastStress.swap(obj.lastStress);
    }

    template<class ScalarType, size_t NumReplica, BaroType Type>
    typename Berendsen<ScalarType, NumReplica, Type>::LatticeMatrix
    Berendsen<ScalarType, NumReplica, Type>::makeDecayMatrix(ScalarType pressPerDOF) const {
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

    template<class ScalarType, size_t NumReplica, BaroType Type>
    typename Berendsen<ScalarType, NumReplica, Type>::LatticeMatrix
    Berendsen<ScalarType, NumReplica, Type>::makeScaleMatrix(
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

    template<class ScalarType, size_t NumReplica, BaroType Type>
    typename Berendsen<ScalarType, NumReplica, Type>::LatticeMatrix
    Berendsen<ScalarType, NumReplica, Type>::makeDeltaLattice(MDType& rpmd, ScalarType deltaT) const {
        const ScalarType pressPerDOF = (tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = makeDecayMatrix(pressPerDOF);
        const auto& lattice = rpmd.getLattice();
        LatticeMatrix result(Dim, Dim, 0);
        auto integrateKernel = [&, deltaT](size_t r, size_t c) {
            using Integrator = ODESolver<ScalarType, 1>;
            using VectorType = typename Integrator::VectorType;
            const VectorType y{result(r, c)};
            result(r, c) = Integrator::rungeKutta4_step(deltaT, 0, y,
                [&lattice, &decayMatrix, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] const VectorType& y) -> VectorType {
                    return {(decayMatrix * lattice).calc(r, c)};
                })[0];
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
}
