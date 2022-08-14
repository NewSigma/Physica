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
    class ReciprocalCell;
    class Poscar;

    class CrystalCell final : public PeriodicCell<Scalar<Float, false>, 3> {
    public:
        using Base = PeriodicCell<Scalar<Float, false>, 3>;
        using ScalarType = Scalar<Float, false>;
        using AtomicArray = Utils::Array<uint16_t>;

        enum class Type : bool {
            Direct,
            Cartesian
        };
    private:
        AtomicArray atomicNumbers;
        Type type;
    public:
        CrystalCell(LatticeMatrix lattice_, PositionMatrix pos_, AtomicArray atomicNumbers_, Type type_);
        CrystalCell(Poscar poscar);
        CrystalCell(const CrystalCell&) = default;
        CrystalCell(CrystalCell&&) noexcept = default;
        ~CrystalCell() = default;
        /* Operators */
        CrystalCell& operator=(CrystalCell cell) noexcept;
        /* Operations */
        void scale(ScalarType factor);
        /* Getters */
        [[nodiscard]] const AtomicArray& getAtomicNumbers() const noexcept { return atomicNumbers; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] size_t getAtomCount() const noexcept { return Base::pos.getRow(); }
        [[nodiscard]] uint16_t getAtomicNumber(size_t ionIndex) const { return atomicNumbers[ionIndex]; }
        [[nodiscard]] ReciprocalCell reciprocal() const noexcept;
        [[nodiscard]] ScalarType getVolume() const noexcept;
        [[nodiscard]] std::unordered_set<uint16_t> getSpecies() const noexcept;
        [[nodiscard]] size_t getElectronCount() const;
        [[nodiscard]] CrystalCell unitToSuper(unsigned int x, unsigned int y, unsigned int z) const;
        /* Helpers */
        void toDirect(PositionMatrix& obj) const;
        void toCartesian(PositionMatrix& obj) const;
        void swap(CrystalCell& cell) noexcept;
    };
}
