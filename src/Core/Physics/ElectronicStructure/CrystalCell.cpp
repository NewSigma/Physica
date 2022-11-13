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
    CrystalCell::CrystalCell(LatticeMatrix lattice_, PositionMatrix pos_, AtomicArray atomicNumbers_, Type type_)
            : Base(std::move(lattice_), std::move(pos_))
            , atomicNumbers(std::move(atomicNumbers_))
            , type(type_) {
        assert(pos.getRow() == atomicNumbers.getLength());
    }

    CrystalCell::CrystalCell(Poscar poscar) : Base(std::move(poscar)), type(poscar.getType()) {}

    CrystalCell& CrystalCell::operator=(CrystalCell cell) noexcept {
        swap(cell);
        return *this;
    }

    void CrystalCell::scale(ScalarType factor) {
        Base::scale_direct(factor);
    }

    void CrystalCell::toDirect() {
        assert(type == Type::Cartesian);
        Base::toDirect(Base::makeInvLattice());
        type = Type::Direct;
    }

    void CrystalCell::toCartesian() {
        assert(type == Type::Direct);
        Base::toCartesian();
        type = Type::Cartesian;
    }

    void CrystalCell::unitToSuper(unsigned int x, unsigned int y, unsigned int z) {
        const size_t numAtom = getAtomCount();
        if (type == Type::Cartesian) {
            toDirect();
            Base::unitToSuper_direct(x, y, z);
            toCartesian();
        }
        else
            Base::unitToSuper_direct(x, y, z);
        
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

    typename CrystalCell::StructureFactorType CrystalCell::makeStructureFactor(ScalarType cutEnergy) const {
        using GridType = KSpaceGrid<ScalarType>;
        using VectorType = Vector<ScalarType, 3>;
        const auto species = getSpecies();
        const auto repCell = Base::reciprocal();
        const auto& lattice = repCell.getLattice();
        auto factorGrids = Utils::Array<GridType>(species.size(), GridType::makeGrid(cutEnergy, lattice));

        size_t i = 0;
        for (uint16_t element : species) {
            GridType& grid = factorGrids[i];
            size_t j = 0;
            GridType::forReducedKInGrid(grid.getDim(), lattice, [this, element, &i, &j, &grid](VectorType k) {
                ComplexType factor = ComplexType::Zero();
                for (size_t ion = 0; ion < getAtomCount(); ++ion) {
                    if (getAtomicNumber(ion) == element) { //Optimize: We can use searching table method
                        const ScalarType phase = k * getPos().row(ion).asVector();
                        ScalarType s, c;
                        sincos(phase, s, c);
                        factor += ComplexType(c, s);
                    }
                }
                grid.asVector()[j] = factor;
                j += 1;
            });
            ++i;
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
        std::swap(type, cell.type);
    }
}
