/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/Container/KSpaceGrid.h"
#include "Physica/Core/IO/Poscar.h"

namespace Physica::Core {
    CrystalCell::CrystalCell(Base base, AtomicArray atomicNumbers_)
            : Base(std::move(base))
            , atomicNumbers(std::move(atomicNumbers_)) {
        assert(pos.getRow() == atomicNumbers.getLength());
    }

    CrystalCell::CrystalCell(Poscar poscar) : Base(std::move(poscar)) {
        assert(poscar.isElementTypeInitialized());
        atomicNumbers.reserve(getNumParticle());
        size_t index = 0;
        for (auto num : poscar.getNumOfEachType()) {
            const uint8_t type = poscar.getElementTypes()[index];
            for (size_t i = 0; i < num; ++i)
                atomicNumbers.append(type);
            index += 1;
        }
    }

    CrystalCell& CrystalCell::operator=(CrystalCell cell) noexcept {
        swap(cell);
        return *this;
    }

    void CrystalCell::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        const size_t numAtom = getAtomCount();
        Base::toSuperCell<ExtendCellOption::DOFMajor>(x, y, z);
        
        const size_t newNumAtom = getAtomCount();
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

    CrystalCell CrystalCell::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const {
        CrystalCell result = *this;
        result.toSuperCell(x, y, z);
        return result;
    }

    typename CrystalCell::StructureFactorType CrystalCell::makeStructureFactor(ScalarType cutEnergy) const {
        using GridType = typename StructureFactorType::ValueType;
        using VectorType = Vector<ScalarType, 3>;
        const auto species = getSpecies();
        const auto repCell = Base::reciprocal();
        const auto& lattice = repCell.getLattice();
        auto factorGrids = Utils::Array<GridType>(species.size(), GridType::makeGrid(cutEnergy, lattice));

        size_t i = 0;
        for (uint16_t element : species) {
            GridType& grid = factorGrids[i];
            size_t j = 0;
            GridType::forReducedKInGrid(grid, lattice, [this, element, &i, &j, &grid](VectorType k) {
                auto factor = ComplexType(0);
                for (size_t ion = 0; ion < getAtomCount(); ++ion) {
                    if (getAtomicNumber(ion) == element) { //Optimize: We can use searching table method
                        const ScalarType phase = k * getPos().row(ion).asVector();
                        ScalarType s, c;
                        sincos(phase, s, c);
                        factor += ComplexType(c, s);
                    }
                }
                grid.flatten()[j] = factor;
                j += 1;
            });
            i += 1;
        }
        return factorGrids;
    }

    std::unordered_set<uint16_t> CrystalCell::getSpecies() const noexcept {
        std::unordered_set<uint16_t> set{};
        for (size_t i = 0; i < atomicNumbers.getLength(); ++i)
            set.emplace(atomicNumbers[i]);
        return set;
    }

    size_t CrystalCell::getElectronCount() const {
        size_t result = 0;
        for (size_t i = 0; i < getAtomCount(); ++i)
            result += getAtomicNumber(i);
        return result;
    }

    void CrystalCell::swap(CrystalCell& cell) noexcept {
        Base::swap(cell);
        atomicNumbers.swap(cell.atomicNumbers);
    }
}
