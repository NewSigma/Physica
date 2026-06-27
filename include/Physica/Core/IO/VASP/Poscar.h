/*
 * Copyright 2021-2026 Weibo He.
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

namespace Physica {
    class Xdatcar;

    template<Scalar T>
    class Poscar final : public PeriodicCell<T, 3> {
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

        using Base = PeriodicCell<T, 3>;
        using LatticeMatrix = Base::LatticeMatrix;
        using Type = Base::Type;
        using ElementTypeArray = Array<uint8_t>;
    private:
        ElementTypeArray elementTypes;
        Array<size_t> numOfEachType;
    public:
        Poscar() = default;
        Poscar(Base base, ElementTypeArray elementTypes_, Array<size_t> numOfEachType_);
        Poscar(CrystalCell<T> cell);
        /* Operators */
        template<Scalar U>
        friend std::ostream& operator<<(std::ostream& os, const Poscar<U>& poscar);
        template<Scalar U>
        friend std::istream& operator>>(std::istream& is, Poscar<U>& poscar);
        /* Operations */
        void standardizeLattice();
        void extendInZ(T factor);
        void toUnitCell(unsigned int x, unsigned int y, unsigned int z);
        void toQECell(std::ostream& os) const;
        void swap(Poscar& __restrict poscar) noexcept;
        /* Getters */
        [[nodiscard]] const Array<uint8_t> getElementTypes() const noexcept { return elementTypes; }
        [[nodiscard]] bool isElementTypeInitialized() const noexcept { return !elementTypes.empty(); }
        [[nodiscard]] const Array<size_t>& getNumOfEachType() const noexcept { return numOfEachType; }
        [[nodiscard]] CrystalSystem getCrystalSystem(double precision) const noexcept;
    private:
        using Base::toUnitCell;
        void readTypesAndNumber(std::istream& is);
        void readAtomPos(std::istream& is);
        [[nodiscard]] size_t sumNumOfEachType() const;
        void extendInZ_direct(T factor);
        /* Static members */
        [[nodiscard]] static uint8_t elementSymbolToNumber(char ch1, char ch2);

        friend class Xdatcar;
    };
}

#include "PoscarImpl.h"
