/*
 * Copyright 2022-2023 WeiBo He.
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

#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim = 3> class MDCell;

    namespace Internal {
        template<class T, unsigned int D>
        class Traits<MDCell<T, D>> {
        public:
            using ScalarType = T;
            constexpr static unsigned int Dim = D;
        };
    }

    template<class ScalarType, unsigned int Dim>
    class MDCell : public PeriodicCell<ScalarType, Dim> {
    public:
        using Base = PeriodicCell<ScalarType, Dim>;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Type;
        using MassVector = Vector<ScalarType>;
    private:
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        explicit MDCell(size_t numParticle);
        template<class OtherScalar>
        MDCell(CrystalCell<OtherScalar> cell);
        template<class OtherScalar>
        MDCell(Poscar<OtherScalar> poscar);
        MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_);
        /* Operations */
        void scale(ScalarType factor);
        void normalize();
        void normalizePos(PositionMatrix& target) const;
        void toDirect(PositionMatrix& target) const { Base::toDirect(target, invLattice); }
        template<ExtendCellOption Option>
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        template<ExtendCellOption Option>
        void toSuperCell(Utils::Array<size_t, 3> size) { toSuperCell<Option>(size[0], size[1], size[2]); }
        template<ExtendCellOption Option>
        MDCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        template<ExtendCellOption Option>
        MDCell makeSuperCell(Utils::Array<size_t, 3> size) const { return makeSuperCell<Option>(size[0], size[1], size[2]); }
        /* Getters */
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * Base::getNumParticle(); }
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] ScalarType getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] const InvLatticeMatrix& getInvLattice() const noexcept { return invLattice; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Type::Cartesian; }
        /* Setters */
        void setLattice(const MDCell& cell);
        void setLattice(LatticeMatrix new_lattice);
    protected:
        void toDirect() { Base::toDirect(invLattice); }
    private:
        using Base::normalize;
        using Base::scale;
        using Base::getType;
        using Base::setLattice;
    };

    template<class ScalarType, unsigned int Dim>
    MDCell<ScalarType, Dim>::MDCell(size_t numParticle) : Base(numParticle), massVec(numParticle) {}

    template<class ScalarType, unsigned int Dim>
    template<class OtherScalar>
    MDCell<ScalarType, Dim>::MDCell(CrystalCell<OtherScalar> cell) {
        if (cell.getType() == Type::Direct)
            cell.toCartesian();
        Base::operator=(Base(cell.getLattice(), cell.getPos(), Type::Cartesian));
        invLattice = Base::makeInvLattice();

        massVec.resize(Base::getNumParticle());
        for (size_t i = 0; i < Base::getNumParticle(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            massVec[i] = PhyConst<AU>::atomMass(atomicNum);
        }
    }

    template<class ScalarType, unsigned int Dim>
    template<class OtherScalar>
    MDCell<ScalarType, Dim>::MDCell(Poscar<OtherScalar> poscar) : MDCell(CrystalCell<OtherScalar>(poscar)) {}

    template<class ScalarType, unsigned int Dim>
    MDCell<ScalarType, Dim>::MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_)
            : Base(std::move(lattice), std::move(pos), Base::Type::Cartesian)
            , massVec(std::move(massVec_)) {
        invLattice = Base::makeInvLattice();
    }

    template<class ScalarType, unsigned int Dim>
    void MDCell<ScalarType, Dim>::scale(ScalarType factor) {
        Base::scale_cartesian(factor);
        invLattice *= Core::reciprocal(factor);
    }

    template<class ScalarType, unsigned int Dim>
    void MDCell<ScalarType, Dim>::normalize() {
        toDirect();
        Base::normalize_direct();
        Base::toCartesian();
    }

    template<class ScalarType, unsigned int Dim>
    void MDCell<ScalarType, Dim>::normalizePos(PositionMatrix& target) const {
        Base::toDirect(target, invLattice);
        for (auto& elem : target) {
            elem -= floor(elem);
            assert(ScalarType(0) <= elem && elem <= ScalarType(1));
        }
        Base::toCartesian(target, Base::getLattice());
    }

    template<class ScalarType, unsigned int Dim>
    template<ExtendCellOption Option>
    void MDCell<ScalarType, Dim>::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        /* Extend mass */ {
            const size_t oldNumParticle = Base::getNumParticle();
            const unsigned int numCell = x * y * z;
            MassVector newMass(oldNumParticle * numCell);        
            if constexpr (Option == ExtendCellOption::DOFMajor) {
                size_t k = 0;
                for (size_t i = 0; i < oldNumParticle; ++i) {
                    for (unsigned int j = 0; j < numCell; ++j) {
                        newMass[k] = massVec[i];
                        k += 1;
                    }
                }
            }
            else {
                size_t k = 0;
                for (unsigned int i = 0; i < numCell; ++i) {
                    for (size_t j = 0; j < oldNumParticle; ++j) {
                        newMass[k] = massVec[j];
                        k += 1;
                    }
                }
            }
            massVec = std::move(newMass);
        }
        toDirect();
        Base::template toSuperCellDirect<Option>(x, y, z);
        Base::toCartesian();
        invLattice = Base::makeInvLattice();
    }

    template<class ScalarType, unsigned int Dim>
    template<ExtendCellOption Option>
    MDCell<ScalarType, Dim>
    MDCell<ScalarType, Dim>::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const {
        MDCell result = *this;
        result.toSuperCell<Option>(x, y, z);
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    void MDCell<ScalarType, Dim>::setLattice(const MDCell& cell) {
        Base::setLattice(cell.getLattice());
        invLattice = cell.invLattice();
    }

    template<class ScalarType, unsigned int Dim>
    void MDCell<ScalarType, Dim>::setLattice(LatticeMatrix new_lattice) {
        Base::setLattice(new_lattice);
        invLattice = Base::makeInvLattice();
    }
}
