/*
 * Copyright 2023-2024 WeiBo He.
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
        inline void npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT);
        using Base::swap;
        /* Getters */
        using Base::getLastStress;
        /* Setters */
        using Base::setTemperature;
    private:
        [[nodiscard]] LatticeMatrix makeDiffuseMatrix(ScalarType pressPerDOF) const;
    };

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::SCRBaro(
            ScalarType compressRate, ScalarType tempT, ScalarType targetP) : Base(compressRate, tempT, targetP) {}

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    template<class MDType, class ForceModel>
    inline void SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::npt_step(
            MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT) {
        const ScalarType pressPerDOF = (Base::tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = Base::makeDecayMatrix(stress, pressPerDOF);
        const auto diffuseMatrix = makeDiffuseMatrix(pressPerDOF);
        Base::scaleRPMD(rpmd, Base::makeDeltaLattice([&, deltaT](size_t r, size_t c) -> ScalarType {
                using Integrator = SRK2<ScalarType, 1>;
                using VectorType = typename Integrator::VectorType;
                [[maybe_unused]] ScalarType unused = 0;
                VectorType sol{ScalarType(0)};
                Integrator::step(deltaT, unused, sol,
                    [&, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] VectorType sol) -> VectorType {
                        return {(decayMatrix * rpmd.getLattice()).calc(r, c)};
                    },
                    [&, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] VectorType sol) -> VectorType {
                        return {(diffuseMatrix * rpmd.getLattice()).calc(r, c)};
                    });
                return sol[0];
            }));
    }

    template<class ScalarType, size_t NumReplica, class RandomPoolType, BaroType Type>
    typename SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::LatticeMatrix
    SCRBaro<ScalarType, NumReplica, RandomPoolType, Type>::makeDiffuseMatrix(ScalarType pressPerDOF) const {
        const ScalarType diffuseFactor = sqrt(ScalarType(2.0 / Dim) * Base::compressRate * pressPerDOF);
        auto& gen = RandomPoolType::getInstance().getGen();
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
}
