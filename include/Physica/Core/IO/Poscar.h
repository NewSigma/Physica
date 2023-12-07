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
#include "Physica/Core/Physics/SolidState/CrystalCell.h"

namespace Physica::Core {
    class Xdatcar;

    template<class ScalarType>
    class Poscar final : public PeriodicCell<ScalarType, 3> {
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

        using Base = PeriodicCell<ScalarType, 3>;
        using LatticeMatrix = typename Base::LatticeMatrix;
        using Type = typename Base::Type;
        using ElementTypeArray = Utils::Array<uint8_t>;
    private:
        using Base::lattice;
        using Base::pos;
        using Base::type;
        ElementTypeArray elementTypes;
        Utils::Array<size_t> numOfEachType;
    public:
        Poscar();
        Poscar(Base base, ElementTypeArray elementTypes_, Utils::Array<size_t> numOfEachType_);
        Poscar(CrystalCell<ScalarType> cell);
        /* Operators */
        template<class AnyScalar>
        friend std::ostream& operator<<(std::ostream& os, const Poscar<AnyScalar>& poscar);
        template<class AnyScalar>
        friend std::istream& operator>>(std::istream& is, Poscar<AnyScalar>& poscar);
        /* Operations */
        void standrizeLattice();
        void extendInZ(ScalarType factor);
        void toUnitCell(unsigned int x, unsigned int y, unsigned int z);
        /* Getters */
        [[nodiscard]] const Utils::Array<uint8_t> getElementTypes() const noexcept { return elementTypes; }
        [[nodiscard]] bool isElementTypeInitialized() const noexcept { return !elementTypes.empty(); }
        [[nodiscard]] const Utils::Array<size_t>& getNumOfEachType() const noexcept { return numOfEachType; }
        [[nodiscard]] CrystalSystem getCrystalSystem(double precision) const noexcept;
        /* Helpers */
        void swap(Poscar& __restrict poscar) noexcept;
    private:
        using Base::toUnitCell;
        void readTypesAndNumber(std::istream& is);
        void readAtomPos(std::istream& is);
        size_t sumNumOfEachType() const;
        void extendInZ_direct(ScalarType factor);
        /* Static members */
        [[nodiscard]] static uint8_t elementSymbolToNumber(char ch1, char ch2);

        friend class Xdatcar;
    };
}

#include "PoscarImpl.h"
