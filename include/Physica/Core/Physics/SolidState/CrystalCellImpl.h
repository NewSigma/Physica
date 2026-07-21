/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/Exception/BadFileFormatException.h"
#include "CrystalCell.h"
#include "VASP/Poscar.h"

namespace Physica {
    template<Scalar T>
    CrystalCell<T>::CrystalCell(Base base, AtomicArray atomicNumbers_)
            : Base(std::move(base))
            , atomicNumbers(std::move(atomicNumbers_)) {
        assert(Base::pos.getRow() == atomicNumbers.getLength());
    }

    template<Scalar T>
    template<Scalar U>
    CrystalCell<T>::CrystalCell(Poscar<U> poscar) : Base(std::move(poscar)) {
        assert(poscar.isElementTypeInitialized());
        atomicNumbers.reserve(Base::getNumParticle());
        size_t index = 0;
        for (auto num : poscar.getNumOfEachType()) {
            const uint8_t type = poscar.getElementTypes()[index];
            for (size_t i = 0; i < num; ++i)
                atomicNumbers.append(type);
            index += 1;
        }
    }

    template<Scalar T>
    void CrystalCell<T>::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        const size_t numAtom = Base::getNumParticle();
        Base::template toSuperCell<ExtendCellOption::AtomMajor>(x, y, z);

        const size_t newNumAtom = Base::getNumParticle();
        AtomicArray new_atomic(newNumAtom);
        size_t index = 0;
        for (size_t atom = 0; atom < numAtom; ++atom) {
            const auto atomic = atomicNumbers[atom];
            for (unsigned int _ = 0; _ < x * y * z; ++_) {
                new_atomic[index] = atomic;
                ++index;
            }
        }
        atomicNumbers.swap(new_atomic);
    }

    template<Scalar T>
    CrystalCell<T> CrystalCell<T>::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const {
        CrystalCell result = *this;
        result.toSuperCell(x, y, z);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    H5Group CrystalCell<T>::read(const H5Loc& loc, const char* name) {
        const auto group = Base::read(loc, name);
        const auto dataset = group.template openDataSet<1>("AtomicNumber");
        const size_t numParticle = dataset.getSize(0);
        if (numParticle != Base::getNumParticle()) [[unlikely]]
            throw BadFileFormatException("[Error]: Bad CrystalCell");
        atomicNumbers.resize(numParticle);
        dataset.read(atomicNumbers.data(), H5::PredType::NATIVE_UINT16);
        return H5Group(group);
    }

    template<Scalar T>
    H5Group CrystalCell<T>::write(H5Loc& loc, const char* name) const {
        auto group = Base::write(loc, name);
        const auto space = H5DataSpace<1>(Base::getNumParticle());
        auto dataset = group.template createDataSet<1>("AtomicNumber", H5::PredType::NATIVE_UINT16, space);
        dataset.write(atomicNumbers.data(), H5::PredType::NATIVE_UINT16);
        return group;
    }
#endif
    template<Scalar T>
    void CrystalCell<T>::swap(CrystalCell& __restrict cell) noexcept {
        assert(this != &cell && "[Error]: Self swap is likely a bug");
        Base::swap(cell);
        atomicNumbers.swap(cell.atomicNumbers);
    }

    template<Scalar T>
    std::unordered_set<uint16_t> CrystalCell<T>::getSpecies() const noexcept {
        std::unordered_set<uint16_t> set{};
        for (size_t i = 0; i < atomicNumbers.getLength(); ++i)
            set.emplace(atomicNumbers[i]);
        return set;
    }

    template<Scalar T>
    size_t CrystalCell<T>::getElectronCount() const {
        const size_t numParticle = Base::getNumParticle();
        size_t result = 0;
        for (size_t i = 0; i < numParticle; ++i)
            result += getAtomicNumber(i);
        return result;
    }
}

