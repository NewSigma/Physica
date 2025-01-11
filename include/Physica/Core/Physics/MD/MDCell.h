/*
 * Copyright 2022-2025 Weibo He.
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

#include <set>
#include <unordered_map>
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    template<Scalar T, unsigned int Dim = 3>
    class MDCell : public PeriodicCell<T, Dim> {
    public:
        using Base = PeriodicCell<T, Dim>;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Type;
        using MassVector = VectorND<T>;
        using ParticleType = uint8_t; // uint8_t should be enough to hold the periodic table
        using MassTypeMap = std::unordered_map<T, ParticleType>;
    private:
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        MDCell() = default;
        explicit MDCell(size_t numParticle);
        template<Scalar U>
        MDCell(CrystalCell<U> cell);
        template<Scalar U>
        MDCell(PeriodicCell<U, Dim> cell, MassVector massVec_);
        template<Scalar U>
        MDCell(Poscar<U> poscar);
        MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_);
        /* Operations */
        void scale(const T& factor);
        void normalize();
        void normalizePos(PositionMatrix& target) const;
        void toDirect(PositionMatrix& target) const { Base::toDirect(target, invLattice); }
        template<ExtendCellOption Option>
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        template<ExtendCellOption Option>
        void toSuperCell(Array<size_t, 3> size) { toSuperCell<Option>(size[0], size[1], size[2]); }
        template<ExtendCellOption Option>
        MDCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        template<ExtendCellOption Option>
        MDCell makeSuperCell(Array<size_t, 3> size) const { return makeSuperCell<Option>(size[0], size[1], size[2]); }

        [[nodiscard]] MassTypeMap makeMassTypeMap() const;
        /* Getters */
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * Base::getNumParticle(); }
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] T getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] const InvLatticeMatrix& getInvLattice() const noexcept { return invLattice; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Type::Cartesian; }
        /* Setters */
        void setLattice(const MDCell& cell);
        void setLattice(LatticeMatrix new_lattice);
        inline void setMass(size_t atom, T value);
    protected:
        void toDirect() { Base::toDirect(invLattice); }
    private:
        using Base::normalize;
        using Base::scale;
        using Base::getType;
        using Base::setLattice;
    };

    template<Scalar T, unsigned int Dim>
    MDCell<T, Dim>::MDCell(size_t numParticle) : Base(numParticle), massVec(numParticle) {}

    template<Scalar T, unsigned int Dim>
    template<Scalar U>
    MDCell<T, Dim>::MDCell(CrystalCell<U> cell) {
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

    template<Scalar T, unsigned int Dim>
    template<Scalar U>
    MDCell<T, Dim>::MDCell(PeriodicCell<U, Dim> cell, MassVector massVec_) : massVec(std::move(massVec_)) {
        if (cell.getType() == Type::Direct)
            cell.toCartesian();
        Base::operator=(Base(cell.getLattice(), cell.getPos(), Type::Cartesian));
        invLattice = Base::makeInvLattice();
    }

    template<Scalar T, unsigned int Dim>
    template<Scalar U>
    MDCell<T, Dim>::MDCell(Poscar<U> poscar) : MDCell(CrystalCell<U>(poscar)) {}

    template<Scalar T, unsigned int Dim>
    MDCell<T, Dim>::MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_)
            : Base(std::move(lattice), std::move(pos), Base::Type::Cartesian)
            , massVec(std::move(massVec_)) {
        invLattice = Base::makeInvLattice();
    }

    template<Scalar T, unsigned int Dim>
    void MDCell<T, Dim>::scale(const T& factor) {
        Base::scale_cartesian(factor);
        invLattice *= Core::reciprocal(factor);
    }

    template<Scalar T, unsigned int Dim>
    void MDCell<T, Dim>::normalize() {
        toDirect();
        Base::normalize_direct();
        Base::toCartesian();
    }

    template<Scalar T, unsigned int Dim>
    void MDCell<T, Dim>::normalizePos(PositionMatrix& target) const {
        Base::toDirect(target, invLattice);
        for (auto& elem : target.asArray()) {
            elem -= floor(elem);
            assert(T(0) <= elem && elem <= T(1));
        }
        Base::toCartesian(target, Base::getLattice());
    }

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    void MDCell<T, Dim>::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        /* Extend mass */ {
            const size_t oldNumParticle = Base::getNumParticle();
            const unsigned int numCell = x * y * z;
            MassVector newMass(oldNumParticle * numCell);        
            if constexpr (Option == ExtendCellOption::AtomMajor) {
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

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    MDCell<T, Dim>
    MDCell<T, Dim>::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const {
        MDCell result = *this;
        result.toSuperCell<Option>(x, y, z);
        return result;
    }

    template<Scalar T, unsigned int Dim>
    MDCell<T, Dim>::MassTypeMap MDCell<T, Dim>::makeMassTypeMap() const {
        ParticleType nextType = 0;
        std::unordered_map<T, ParticleType> massToType{};
        std::set<T> massSet{};
        for (auto mass : massVec)
            massSet.insert(mass);

        for (T mass : massSet) {
            massToType.insert(std::make_pair(mass, nextType));
            nextType += 1;
        }
        return massToType;
    }

    template<Scalar T, unsigned int Dim>
    void MDCell<T, Dim>::setLattice(const MDCell& cell) {
        Base::setLattice(cell.getLattice());
        invLattice = cell.invLattice();
    }

    template<Scalar T, unsigned int Dim>
    void MDCell<T, Dim>::setLattice(LatticeMatrix new_lattice) {
        Base::setLattice(new_lattice);
        invLattice = Base::makeInvLattice();
    }

    template<Scalar T, unsigned int Dim>
    inline void MDCell<T, Dim>::setMass(size_t atom, T value) {
        assert(atom < Base::getNumParticle() && "[Error]: Index out of range");
        assert(value.isPositive() && "[Error]: Mass must be positive");
        massVec[atom] = value;
    }
}

namespace Physica {
    template<Scalar T, unsigned int D>
    class Traits<Core::MDCell<T, D>> {
    public:
        using ScalarType = T;
        constexpr static unsigned int Dim = D;
    };
}
