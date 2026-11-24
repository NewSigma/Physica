/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Physics/MD/Barostat/Berendsen.h"
#include "FireModel.h"

namespace Physica {
    /**
     * \class CFireModel is Cell-FIRE that incorporates cell degree of freedom into \class FireModel
     */
    template<Scalar T, unsigned int Dim, BaroType Type>
    class CFireModel : private FireModel<T, Dim> {
        using This = CFireModel<T, Dim, Type>;
        using Base = FireModel<T, Dim>;
        using typename Base::MDType;
        using LatticeMatrix = MDType::LatticeMatrix;
        using BerendsenType = Berendsen<T, 1, Type>;

        BerendsenType baro;
        LatticeMatrix lattMomentum;
        LatticeMatrix deltaLattice;
        T lattMass;
        T normLattF;
    public:
        CFireModel(Base fire, T targetP, T lattMass_);
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
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getTimeStep;
        using Base::getForceNorm;
        [[nodiscard]] const LatticeMatrix& getLattMomentum() const noexcept { return lattMomentum; }
        [[nodiscard]] const LatticeMatrix& getLastStress() const noexcept { return baro.getLastStress(); }
        [[nodiscard]] T getLattMass() const noexcept { return lattMass; }
        [[nodiscard]] T getNormLattF() const noexcept { return normLattF; }
    };

    template<Scalar T, unsigned int Dim, BaroType Type>
    CFireModel<T, Dim, Type>::CFireModel(Base fire, T targetP, T lattMass_)
            : Base(std::move(fire))
            , baro(reciprocal(lattMass_), 0, targetP)
            , lattMass(lattMass_) {
        lattMomentum.zeros();
    }

    template<Scalar T, unsigned int Dim, BaroType Type>
    void CFireModel<T, Dim, Type>::nve_step(MDType& rpmd) const {
        const LatticeMatrix deltaLattice = lattMomentum * (rpmd.getTimeStep() / lattMass);
        BerendsenType::scaleRPMD(rpmd, deltaLattice);
    }

    template<Scalar T, unsigned int Dim, BaroType Type>
    void CFireModel<T, Dim, Type>::forceStep(MDType& rpmd, LatticeMatrix stress) {
        const LatticeMatrix decayMatrix = baro.makeDecayMatrix(stress, 0);
        deltaLattice = BerendsenType::makeDeltaLattice([this, &rpmd, &decayMatrix](size_t r, size_t c) -> T {
                return (decayMatrix * rpmd.getLattice()).calc(r, c) * getTimeStep();
            });
        lattMomentum += deltaLattice;
    }

    template<Scalar T, unsigned int Dim, BaroType Type>
    void CFireModel<T, Dim, Type>::mixingStep(MDType& rpmd) {
        Base::mixingStep(rpmd);
        const T normLattP = lattMomentum.flatten().norm();
        normLattF = deltaLattice.flatten().norm();
        const T mixAlpha = Base::getMixAlpha();
        lattMomentum = (T(1) - mixAlpha) * lattMomentum + (mixAlpha * normLattP / normLattF) * deltaLattice;
    }

    template<Scalar T, unsigned int Dim, BaroType Type>
    void CFireModel<T, Dim, Type>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        baro.swap(obj.baro);
        lattMomentum.swap(obj.lattMomentum);
        lattMass.swap(obj.lattMass);
        normLattF.swap(obj.normLattF);
    }
}
