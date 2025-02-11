/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica {
    template<Scalar T> class Poscar;

    template<Scalar T>
    class CrystalCell final : public PeriodicCell<T, 3> {
    public:
        constexpr static unsigned int Dim = 3;
        using Base = PeriodicCell<T, Dim>;
        using ComplexType = Complex<T>;
        using AtomicArray = Array<uint16_t>;
        using typename Base::Type;
    private:
        AtomicArray atomicNumbers;
    public:
        CrystalCell() = default;
        CrystalCell(Base base, AtomicArray atomicNumbers_);
        template<Scalar U>
        CrystalCell(Poscar<U> poscar);
        CrystalCell(const CrystalCell&) = default;
        CrystalCell(CrystalCell&&) noexcept = default;
        ~CrystalCell() = default;
        /* Operators */
        CrystalCell& operator=(CrystalCell cell) noexcept;
        /* Operations */
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        [[nodiscard]] CrystalCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(CrystalCell& __restrict cell) noexcept;
        /* Getters */
        using Base::getType;
        [[nodiscard]] const AtomicArray& getAtomicNumbers() const noexcept { return atomicNumbers; }
        [[nodiscard]] uint16_t getAtomicNumber(size_t ionIndex) const { return atomicNumbers[ionIndex]; }
        [[nodiscard]] std::unordered_set<uint16_t> getSpecies() const noexcept;
        [[nodiscard]] size_t getElectronCount() const;
    private:
        using Base::toSuperCell;
    };
}

#include "CrystalCellImpl.h"
