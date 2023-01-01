/*
 * Copyright 2021 WeiBo He.
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

#include <iosfwd>
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"

namespace Physica::Core {
    class Xdatcar;
    class CrystalCell;

    class Poscar final : public PeriodicCell<typename CrystalCell::ScalarType, 3> {
    public:
        enum CrystalSystem : char {
            Triclinic,
            Monoclinic,
            Orthohombic,
            Tetragonal,
            Hexagonal,
            Rhombohedral,
            Cubic
        };

        using Base = PeriodicCell<typename CrystalCell::ScalarType, 3>;
        using ScalarType = typename CrystalCell::ScalarType;
        using Type = typename CrystalCell::Type;
    private:
        Utils::Array<size_t> numOfEachType;
        Type type;
    public:
        Poscar();
        Poscar(LatticeMatrix lattice_, PositionMatrix pos_, Utils::Array<size_t> numOfEachType_, Type type_);
        Poscar(CrystalCell cell);
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Poscar& poscar);
        friend std::istream& operator>>(std::istream& is, Poscar& poscar);
        /* Operations */
        void scale(ScalarType factor);
        void standrizeLattice();
        void extendInZ(ScalarType factor);
        void superToUnit(unsigned int x, unsigned int y, unsigned int z);
        /* Getters */
        [[nodiscard]] const Utils::Array<size_t>& getNumOfEachType() const noexcept { return numOfEachType; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] CrystalSystem getCrystalSystem(double precision) const noexcept;
        [[nodiscard]] size_t getAtomCount() const noexcept { return pos.getRow(); }
        /* Helpers */
        void swap(Poscar& poscar) noexcept;
    private:
        void readNumOfEachType(std::istream& is);
        void readAtomPos(std::istream& is);
        size_t sumNumOfEachType() const;

        friend class Xdatcar;
    };
}

namespace std {
    template<>
    inline void swap<Physica::Core::Poscar>(
            Physica::Core::Poscar& car1, Physica::Core::Poscar& car2) noexcept {
        car1.swap(car2);
    }
}
