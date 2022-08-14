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
        assert(factor.isPositive());
        lattice *= factor;
    }

    typename CrystalCell::ScalarType CrystalCell::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2).asVector());
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

    CrystalCell CrystalCell::unitToSuper(unsigned int x, unsigned int y, unsigned int z) const {
        assert(x > 0 && y > 0 && z > 0);
        const size_t numAtom = x * y * z * getAtomCount();
        LatticeMatrix new_lattice{};
        {
            auto row1 = new_lattice.row(0);
            row1 = lattice.row(0).asVector() * ScalarType(x);
            auto row2 = new_lattice.row(1);
            row2 = lattice.row(1).asVector() * ScalarType(y);
            auto row3 = new_lattice.row(2);
            row3 = lattice.row(2).asVector() * ScalarType(z);
        }
        PositionMatrix new_pos(numAtom, 3);
        AtomicArray new_atomic(numAtom);
        size_t index = 0;
        for (size_t atom = 0; atom < getAtomCount(); ++atom) {
            const auto atomic = atomicNumbers[atom];
            for (unsigned int z_ = 0; z_ < z; ++z_) {
                for (unsigned int y_ = 0; y_ < y; ++y_) {
                    for (unsigned int x_ = 0; x_ < x; ++x_) {
                        new_pos(index, 0) = (pos(atom, 0) + ScalarType(x_)) / ScalarType(x);
                        new_pos(index, 1) = (pos(atom, 1) + ScalarType(y_)) / ScalarType(y);
                        new_pos(index, 2) = (pos(atom, 2) + ScalarType(z_)) / ScalarType(z);
                        new_atomic[index] = atomic;
                        ++index;
                    }
                }
            }
        }
        return CrystalCell(std::move(new_lattice), std::move(new_pos), std::move(new_atomic), type);
    }

    void CrystalCell::toDirect(PositionMatrix& obj) const {
        assert(type == Type::Cartesian);
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, 3, 3>;
        const MatrixType inv = lattice.inverse();
        obj *= inv;
    }

    void CrystalCell::toCartesian(PositionMatrix& obj) const {
        assert(type == Type::Direct);
        obj *= lattice;
    }

    void CrystalCell::swap(CrystalCell& cell) noexcept {
        Base::swap(cell);
        atomicNumbers.swap(cell.atomicNumbers);
        std::swap(type, cell.type);
    }
}
