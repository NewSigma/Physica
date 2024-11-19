/*
 * Copyright 2023-2024 Weibo He.
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
    /**
     * Stochastic cell rescaling(SCR) barostat as introduced in [1], its anisotropic version as introduced in [2]
     * 
     * Reference:
     * [1] J. Chem. Phys. 153, 114107 (2020); https://doi.org/10.1063/5.0020514
     * [2] https://doi.org/10.48550/arXiv.2111.06402
     */
    template<Scalar T, size_t NumReplica, class RandomType, BaroType Type>
    class SCRBaro : private Berendsen<T, NumReplica, Type> {
        using Base = Berendsen<T, NumReplica, Type>;
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using Base::Dim;
    public:
        SCRBaro() = default;
        SCRBaro(T compressRate, T tempT, T targetP);
        SCRBaro(const SCRBaro&) = default;
        SCRBaro(SCRBaro&&) noexcept = default;
        ~SCRBaro() = default;
        /* Operators */
        SCRBaro& operator=(SCRBaro obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class MDType, class ForceModel>
        inline void npt_step(MDType& rpmd, const LatticeMatrix& stress, T deltaT);
        using Base::swap;
        /* Getters */
        using Base::getLastStress;
        /* Setters */
        using Base::setTemperature;
    private:
        [[nodiscard]] LatticeMatrix makeDiffuseMatrix(T pressPerDOF) const;
    };

    template<Scalar T, size_t NumReplica, class RandomType, BaroType Type>
    SCRBaro<T, NumReplica, RandomType, Type>::SCRBaro(
            T compressRate, T tempT, T targetP) : Base(compressRate, tempT, targetP) {}

    template<Scalar T, size_t NumReplica, class RandomType, BaroType Type>
    template<class MDType, class ForceModel>
    inline void SCRBaro<T, NumReplica, RandomType, Type>::npt_step(
            MDType& rpmd, const LatticeMatrix& stress, T deltaT) {
        const T pressPerDOF = (Base::tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = Base::makeDecayMatrix(stress, pressPerDOF);
        const auto diffuseMatrix = makeDiffuseMatrix(pressPerDOF);
        Base::scaleRPMD(rpmd, Base::makeDeltaLattice([&, deltaT](size_t r, size_t c) -> T {
                using Integrator = SRK2<T, 1>;
                using VectorType = typename Integrator::VectorType;
                [[maybe_unused]] T unused = 0;
                VectorType sol{T(0)};
                Integrator::step(deltaT, unused, sol,
                    [&, r, c]([[maybe_unused]] T x, [[maybe_unused]] VectorType sol) -> VectorType {
                        return {(decayMatrix * rpmd.getLattice()).calc(r, c)};
                    },
                    [&, r, c]([[maybe_unused]] T x, [[maybe_unused]] VectorType sol) -> VectorType {
                        return {(diffuseMatrix * rpmd.getLattice()).calc(r, c)};
                    });
                if (abs(sol[0]) > abs(rpmd.getLattice().row(r)).max()) [[unlikely]]
                    throw std::runtime_error("[Error]: Unexpected large delta lattice");
                return sol[0];
            }));
    }

    template<Scalar T, size_t NumReplica, class RandomType, BaroType Type>
    typename SCRBaro<T, NumReplica, RandomType, Type>::LatticeMatrix
    SCRBaro<T, NumReplica, RandomType, Type>::makeDiffuseMatrix(T pressPerDOF) const {
        const T diffuseFactor = sqrt(T(2.0 / Dim) * Base::compressRate * pressPerDOF);
        auto& gen = RandomType::getInstance().getGen();
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (Type == BaroType::Anisotropic) {
            const auto rand = LatticeMatrix::random_normal(Dim, Dim, gen);
            result = (rand + rand.transpose()) * (diffuseFactor * T(0.5));
        }
        else if constexpr (Type == BaroType::XY) {
            result(0, 0) = T::random_normal(gen);
            result(1, 1) = T::random_normal(gen);
            result(0, 1) = result(1, 0) = T::random_normal(gen) * T(M_SQRT1_2);
            auto corner = result.topLeftCorner(2);
            corner *= diffuseFactor;
        }
        else if constexpr (Type == BaroType::Z) {
            result(2, 2) = T::random_normal(gen) * diffuseFactor;
        }
        else {
            constexpr bool Unreachable = Type == BaroType::Anisotropic;
            static_assert(Unreachable, "[Error]: Not implemented");
        }
        return result;
    }
}

namespace Physica {
    template<Scalar T, size_t NumReplica, class RandomType, BaroType Type>
    class Traits<Core::SCRBaro<T, NumReplica, RandomType, Type>> {
    public:
        constexpr static unsigned int Order = 1;
    };
}
