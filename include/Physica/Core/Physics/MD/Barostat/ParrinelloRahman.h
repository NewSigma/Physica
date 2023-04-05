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

namespace Physica::Core {
    /**
     * Reference:
     * [1] M. Parrinello and A. Rahman, J. Appl. Phys. 52, 7182 (1981); doi: 10.1063/1.328693
     */
    template<class ScalarType, class PosScalarType>
    class ParrinelloRahman {
        using MDType = RPMD<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDType::MDCellType::LatticeMatrix;

        ScalarType latticeMass;
        LatticeMatrix latticeMomentum;
        LatticeMatrix targetStress;
        
        LatticeMatrix buffer;
    public:
        ParrinelloRahman() = default;
        ParrinelloRahman(ScalarType latticeMass_, LatticeMatrix targetStress_);
        ParrinelloRahman(const ParrinelloRahman&) = default;
        ParrinelloRahman(ParrinelloRahman&&) noexcept = default;
        ~ParrinelloRahman() = default;
        /* Operators */
        ParrinelloRahman& operator=(ParrinelloRahman obj) noexcept;
        /* Operations */
        void forceStep(MDType& md, ScalarType deltaT);
        void swap(ParrinelloRahman& obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getLatticeMass() const noexcept { return latticeMass; }
        [[nodiscard]] const LatticeMatrix& getLatticeMomentum() const noexcept { return latticeMomentum; }
    };

    template<class ScalarType, class PosScalarType>
    ParrinelloRahman<ScalarType, PosScalarType>::ParrinelloRahman(ScalarType latticeMass_, LatticeMatrix targetStress_)
            : latticeMass(latticeMass_)
            , latticeMomentum(3, 3, 0)
            , targetStress(std::move(targetStress_)) {
        buffer = targetStress.transpose();
        targetStress = (targetStress + buffer) * ScalarType(0.5);
    }

    template<class ScalarType, class PosScalarType>
    ParrinelloRahman<ScalarType, PosScalarType>&
    ParrinelloRahman<ScalarType, PosScalarType>::operator=(ParrinelloRahman obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    void ParrinelloRahman<ScalarType, PosScalarType>::forceStep(MDType& md, ScalarType deltaT) {
        assert(md.getNumReplica() == 1);
        constexpr unsigned int Dim = md.getDim();
        /* Make lattice correction matrix */ {
            const ScalarType factor = reciprocal(ScalarType(2 * M_PI) * latticeMass);
            buffer = md.getInvLattice().transpose() * latticeMomentum * factor;
        }
        const size_t numReplica = md.getNumReplica();
        const size_t dof = md.getDOF();
        const size_t numParticle = md.getNumParticle();
        auto& phaseMatrix = md.getPhaseMatrix();
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto col = phaseMatrix.col(replica);
            auto momentum = col.head(dof);
            auto force = md.getForce().col(replica);
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                auto p = momentum.segment(from, to);
                const Vector<ScalarType, 3> deltaP = (force.segment(from, to) - buffer * p)  * deltaT;
                p += deltaP;
            }
        }

        const ScalarType factor = deltaT / ScalarType(2 * M_PI);
        const LatticeMatrix stress = md.makeStress();
        buffer = (stress + stress.transpose()) * ScalarType(0.5);
        latticeMomentum += factor * (md.getInvLattice().transpose() * (buffer - targetStress)).compute();
    }

    template<class ScalarType, class PosScalarType>
    void ParrinelloRahman<ScalarType, PosScalarType>::swap(ParrinelloRahman& obj) noexcept {
        latticeMass.swap(obj.latticeMass);
        latticeMomentum.swap(obj.latticeMomentum);
        targetStress.swap(obj.targetStress);
        buffer.swap(obj.buffer);
    }
}
