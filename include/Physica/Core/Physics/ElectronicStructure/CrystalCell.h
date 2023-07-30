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

#include <unordered_set>
#include "Physica/Core/Physics/PeriodicCell.h"

namespace Physica::Core {
    class Poscar;
    template<class T> class KSpaceGrid;

    class CrystalCell final : public PeriodicCell<Scalar<Float>, 3> {
    public:
        using Base = PeriodicCell<Scalar<Float>, 3>;
        using ScalarType = Scalar<Float>;
        using ComplexType = ComplexScalar<ScalarType>;
        using AtomicArray = Utils::Array<uint16_t>;
        using StructureFactorType = Utils::Array<KSpaceGrid<ComplexType>>;
    private:
        AtomicArray atomicNumbers;
    public:
        CrystalCell() = default;
        CrystalCell(Base base, AtomicArray atomicNumbers_);
        CrystalCell(Poscar poscar);
        CrystalCell(const CrystalCell&) = default;
        CrystalCell(CrystalCell&&) noexcept = default;
        ~CrystalCell() = default;
        /* Operators */
        CrystalCell& operator=(CrystalCell cell) noexcept;
        /* Operations */
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        [[nodiscard]] CrystalCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        [[nodiscard]] StructureFactorType makeStructureFactor(ScalarType cutEnergy) const;
        /* Getters */
        [[nodiscard]] const AtomicArray& getAtomicNumbers() const noexcept { return atomicNumbers; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] size_t getAtomCount() const noexcept { return Base::pos.getRow(); }
        [[nodiscard]] uint16_t getAtomicNumber(size_t ionIndex) const { return atomicNumbers[ionIndex]; }
        [[nodiscard]] std::unordered_set<uint16_t> getSpecies() const noexcept;
        [[nodiscard]] size_t getElectronCount() const;
        /* Helpers */
        void swap(CrystalCell& cell) noexcept;
    private:
        using Base::toSuperCell;
    };
}
