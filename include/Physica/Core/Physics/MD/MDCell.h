/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/PeriodicCell.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType>
    class MDCell : public PeriodicCell<PosScalarType, 3> {
    public:
        constexpr static unsigned int Dim = 3;
        using Base = PeriodicCell<PosScalarType, 3>;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Type;
        using MassVector = Vector<ScalarType>;
    private:
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        MDCell(CrystalCell cell);
        MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_);
        /* Operations */
        void scale(PosScalarType factor);
        void normalize();
        void normalizePos(PositionMatrix& target) const;
        void toDirect(PositionMatrix& target) const { Base::toDirect(target, invLattice); }
        void unitToSuper(unsigned int x, unsigned int y, unsigned int z);
        /* Getters */
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] ScalarType getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] const InvLatticeMatrix& getInvLattice() const noexcept { return invLattice; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Type::Cartesian; }
    protected:
        void toDirect() { Base::toDirect(invLattice); }
    private:
        using Base::normalize;
        using Base::scale;
        using Base::getType;
    };

    template<class ScalarType, class PosScalarType>
    MDCell<ScalarType, PosScalarType>::MDCell(CrystalCell cell) {
        if (cell.getType() == Type::Direct)
            cell.toCartesian();
        Base::operator=(Base(cell.getLattice(), cell.getPos(), Type::Cartesian));
        invLattice = Base::makeInvLattice();

        massVec.resize(Base::getNumParticle());
        for (size_t i = 0; i < Base::getNumParticle(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            massVec[i] = PhyConst<AU>::atomMass(atomicNum);
        }
        normalize();
    }

    template<class ScalarType, class PosScalarType>
    MDCell<ScalarType, PosScalarType>::MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_)
            : Base(std::move(lattice), std::move(pos), Base::Type::Cartesian)
            , massVec(std::move(massVec_)) {
        invLattice = Base::makeInvLattice();
        normalize();
    }

    template<class ScalarType, class PosScalarType>
    void MDCell<ScalarType, PosScalarType>::scale(PosScalarType factor) {
        Base::scale_cartesian(factor);
        invLattice *= Core::reciprocal(factor);
    }

    template<class ScalarType, class PosScalarType>
    void MDCell<ScalarType, PosScalarType>::normalize() {
        Base::toDirect(invLattice);
        Base::normalize_direct();
        Base::toCartesian();
    }

    template<class ScalarType, class PosScalarType>
    void MDCell<ScalarType, PosScalarType>::normalizePos(PositionMatrix& target) const {
        Base::toDirect(target, invLattice);
        for (auto& elem : target) {
            elem -= floor(elem);
            assert(PosScalarType::Zero() <= elem && elem <= PosScalarType::One());
        }
        Base::toCartesian(target, Base::getLattice());
    }

    template<class ScalarType, class PosScalarType>
    void MDCell<ScalarType, PosScalarType>::unitToSuper(unsigned int x, unsigned int y, unsigned int z) {
        toDirect();
        Base::unitToSuper_direct(x, y, z);
        Base::toCartesian();
    }
}
