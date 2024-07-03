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

#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "BaroType.h"

namespace Physica::Core {
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
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using InvLatticeMatrix = typename MDCellType::InvLatticeMatrix;
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
        [[nodiscard]] LatticeMatrix makeDecayMatrix(const LatticeMatrix& stress, ScalarType pressPerDOF);
        template<class MDType, class ForceModel>
        void npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT);
        void swap(Berendsen& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLastStress() const noexcept { return lastStress; }
        /* Setters */
        void setTemperature(ScalarType tempT_) { tempT = std::move(tempT_); }
        /* Static members */
        [[nodiscard]] static LatticeMatrix makeScaleMatrix(const InvLatticeMatrix& invLatt, const LatticeMatrix& deltaLattice);
        template<class MDType>
        static void scaleRPMD(MDType& rpmd, const LatticeMatrix& deltaLattice);
        template<class Integrator>
        [[nodiscard]] static LatticeMatrix makeDeltaLattice(Integrator kernel);
    };

    template<class ScalarType, size_t NumReplica, BaroType Type>
    Berendsen<ScalarType, NumReplica, Type>::Berendsen(
            ScalarType compressRate_, ScalarType tempT_, ScalarType targetP_)
            : compressRate(compressRate_)
            , tempT(tempT_)
            , targetP(targetP_)
            , lastStress() {
        assert(compressRate.isPositive() && "[Error]: Compress rate must be positive");
        assert(!tempT.isNegative() && "[Error]: Temperature must not be negative");
    }

    template<class ScalarType, size_t NumReplica, BaroType Type>
    typename Berendsen<ScalarType, NumReplica, Type>::LatticeMatrix
    Berendsen<ScalarType, NumReplica, Type>::makeDecayMatrix(const LatticeMatrix& stress, ScalarType pressPerDOF) {
        lastStress = stress;
        LatticeMatrix result{};
        result = ScalarType(0);
        const ScalarType centerTargetP = targetP - pressPerDOF;
        if constexpr (Type == BaroType::Anisotropic) {
            result = lastStress;
            for (size_t i = 0; i < Dim; ++i)
                result(i, i) -= centerTargetP;
        }
        else if constexpr (Type == BaroType::XY) {
            auto corner = result.topLeftCorner(2);
            corner = lastStress.topLeftCorner(2);
            for (size_t i = 0; i < 2; ++i)
                result(i, i) -= centerTargetP;
        }
        else if constexpr (Type == BaroType::Z) {
            result(2, 2) = lastStress(2, 2) - centerTargetP;
        }
        else {
            constexpr bool Unreachable = Type == BaroType::Anisotropic;
            static_assert(Unreachable, "[Error]: Not implemented");
        }
        result *= compressRate / ScalarType(Dim);
        return result;
    }

    template<class ScalarType, size_t NumReplica, BaroType Type>
    template<class MDType, class ForceModel>
    void Berendsen<ScalarType, NumReplica, Type>::npt_step(MDType& rpmd, const LatticeMatrix& stress, ScalarType deltaT) {
        const ScalarType pressPerDOF = (tempT * PhyConst<AU>::boltzmannK) / rpmd.getVolume();
        const auto decayMatrix = makeDecayMatrix(stress, pressPerDOF);
        scaleRPMD(rpmd, makeDeltaLattice([&, deltaT](size_t r, size_t c) -> ScalarType {
                using Integrator = ODESolver<ScalarType, 1>;
                using VectorType = typename Integrator::VectorType;
                const VectorType y{ScalarType(0)};
                const ScalarType deltaElem = Integrator::rungeKutta4_step(deltaT, 0, y,
                    [&, r, c]([[maybe_unused]] ScalarType x, [[maybe_unused]] const VectorType& y) -> VectorType {
                        return {(decayMatrix * rpmd.getLattice()).calc(r, c)};
                    })[0];
                if (abs(deltaElem) > abs(rpmd.getLattice()(r, c))) [[unlikely]]
                    throw std::runtime_error("[Error]: Unexpected large delta lattice");
                return deltaElem;
            }));
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
    Berendsen<ScalarType, NumReplica, Type>::makeScaleMatrix(const InvLatticeMatrix& invLatt, const LatticeMatrix& deltaLattice) {
        LatticeMatrix result;
        if constexpr (Type == BaroType::Anisotropic)
            result = deltaLattice * invLatt;
        else {
            result = deltaLattice * invLatt;
            auto col = result.col(2);
            col = ScalarType(0);
        }
        return result;
    }

    template<class ScalarType, size_t NumReplica, BaroType Type>
    template<class MDType>
    void Berendsen<ScalarType, NumReplica, Type>::scaleRPMD(MDType& rpmd, const LatticeMatrix& deltaLattice) {
        const size_t numReplica = rpmd.getNumReplica();
        const size_t numParticle = rpmd.getNumParticle();
        const size_t dof = rpmd.getDOF();
        auto& phase = rpmd.getPhaseMatrix();
        const LatticeMatrix scaleMatrix = makeScaleMatrix(rpmd.getInvLattice(), deltaLattice);
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
    template<class Integrator>
    typename Berendsen<ScalarType, NumReplica, Type>::LatticeMatrix
    Berendsen<ScalarType, NumReplica, Type>::makeDeltaLattice(Integrator kernel) {
        using ResultType = typename std::invoke_result<Integrator, size_t, size_t>::type;
        static_assert(std::is_same<ScalarType, ResultType>::value, "[Error]: Invalid integrator");
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (Type == BaroType::Anisotropic) {
            for (size_t r = 0; r < Dim; ++r)
                for (size_t c = 0; c < Dim; ++c)
                    result(r, c) = kernel(r, c);
        }
        else if constexpr (Type == BaroType::XY) {
            for (size_t r = 0; r < 2; ++r)
                for (size_t c = 0; c < 2; ++c)
                    result(r, c) = kernel(r, c);
        }
        else if constexpr (Type == BaroType::Z) {
            result(2, 2) = kernel(2, 2);
        }
        else {
            constexpr bool Unreachable = Type == BaroType::Anisotropic;
            static_assert(Unreachable, "[Error]: Not implemented");
        }
        return result;
    }
}

namespace Physica {
    template<class ScalarType, size_t NumReplica, BaroType Type>
    class Traits<Core::Berendsen<ScalarType, NumReplica, Type>> {
    public:
        constexpr static unsigned int Order = 1;
    };
}
