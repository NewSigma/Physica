/*
 * Copyright 2023 Weibo He.
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

#include "FireModel.h"
#include "Physica/Core/Physics/MD/Barostat/Berendsen.h"

namespace Physica::Core {
    /**
     * \class CFireModel is Cell-FIRE that incorporates cell degree of freedom into \class FireModel
     */
    template<class ScalarType, unsigned int Dim, BaroType Type>
    class CFireModel : private FireModel<ScalarType, Dim> {
        using Base = FireModel<ScalarType, Dim>;
        using typename Base::MDType;
        using LatticeMatrix = typename MDType::LatticeMatrix;
        using BerendsenType = Berendsen<ScalarType, 1, Type>;

        BerendsenType baro;
        LatticeMatrix lattMomentum;
        LatticeMatrix deltaLattice;
        ScalarType lattMass;
        ScalarType normLattF;
    public:
        CFireModel(Base fire, ScalarType targetP, ScalarType lattMass_);
        CFireModel(const CFireModel&) = default;
        CFireModel(CFireModel&&) noexcept = default;
        ~CFireModel() = default;
        /* Operators */
        CFireModel& operator=(CFireModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::paramStep;
        void nve_step(MDType& rpmd) const;
        void forceStep(MDType& rpmd, LatticeMatrix stress);
        void mixingStep(MDType& rpmd);
        void swap(CFireModel& __restrict obj) noexcept;
        /* Getters */
        using Base::getTimeStep;
        using Base::getForceNorm;
        [[nodiscard]] const LatticeMatrix& getLattMomentum() const noexcept { return lattMomentum; }
        [[nodiscard]] const LatticeMatrix& getLastStress() const noexcept { return baro.getLastStress(); }
        [[nodiscard]] ScalarType getLattMass() const noexcept { return lattMass; }
        [[nodiscard]] ScalarType getNormLattF() const noexcept { return normLattF; }
    };

    template<class ScalarType, unsigned int Dim, BaroType Type>
    CFireModel<ScalarType, Dim, Type>::CFireModel(Base fire, ScalarType targetP, ScalarType lattMass_)
            : Base(std::move(fire))
            , baro(reciprocal(lattMass_), 0, targetP)
            , lattMass(lattMass_) {
        lattMomentum = ScalarType(0);
    }

    template<class ScalarType, unsigned int Dim, BaroType Type>
    void CFireModel<ScalarType, Dim, Type>::nve_step(MDType& rpmd) const {
        const LatticeMatrix deltaLattice = lattMomentum * (rpmd.getTimeStep() / lattMass);
        BerendsenType::scaleRPMD(rpmd, deltaLattice);
    }

    template<class ScalarType, unsigned int Dim, BaroType Type>
    void CFireModel<ScalarType, Dim, Type>::forceStep(MDType& rpmd, LatticeMatrix stress) {
        const LatticeMatrix decayMatrix = baro.makeDecayMatrix(stress, 0);
        deltaLattice = BerendsenType::makeDeltaLattice([this, &rpmd, &decayMatrix](size_t r, size_t c) -> ScalarType {
                return (decayMatrix * rpmd.getLattice()).calc(r, c) * getTimeStep();
            });
        lattMomentum += deltaLattice;
    }

    template<class ScalarType, unsigned int Dim, BaroType Type>
    void CFireModel<ScalarType, Dim, Type>::mixingStep(MDType& rpmd) {
        Base::mixingStep(rpmd);
        const ScalarType normLattP = lattMomentum.flatten().norm();
        normLattF = deltaLattice.flatten().norm();
        const ScalarType mixAlpha = Base::getMixAlpha();
        lattMomentum = (ScalarType(1) - mixAlpha) * lattMomentum + (mixAlpha * normLattP / normLattF) * deltaLattice;
    }

    template<class ScalarType, unsigned int Dim, BaroType Type>
    void CFireModel<ScalarType, Dim, Type>::swap(CFireModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        baro.swap(obj.baro);
        lattMomentum.swap(obj.lattMomentum);
        lattMass.swap(obj.lattMass);
        normLattF.swap(obj.normLattF);
    }
}
