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
#include <iostream>
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

namespace Physica::Core {
    Poscar::Poscar() : lattice(LatticeMatrix::unitMatrix(3))
                     , pos()
                     , numOfEachType()
                     , type(Direct) {}

    std::ostream& operator<<(std::ostream& os, const Poscar& poscar) {
        os << '\n';
        os << 1.0 << '\n';
        os << poscar.lattice << '\n';
        for (size_t i = 0; i < poscar.numOfEachType.getLength(); ++i)
            os << '\t' << poscar.numOfEachType[i];
        os << '\n';
        os << ((poscar.type == Poscar::Direct) ? "Direct\n" : "Cartesian\n");
        os << poscar.pos;
        os << '\n';
        return os;
    }

    std::istream& operator>>(std::istream& is, Poscar& poscar) {
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        is >> poscar.lattice(0, 0) >> poscar.lattice(0, 1) >> poscar.lattice(0, 2);
        is >> poscar.lattice(1, 0) >> poscar.lattice(1, 1) >> poscar.lattice(1, 2);
        is >> poscar.lattice(2, 0) >> poscar.lattice(2, 1) >> poscar.lattice(2, 2);
        /* Test the format of poscar */ {
            int temp;
            is >> temp;
            if (is.fail()) {
                is.clear();
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        /* Read numOfEachType */ {
            poscar.numOfEachType.reserve(8);
            do {
                size_t count;
                is >> count;
                if (is.good())
                    poscar.numOfEachType << count;
                else
                    break;
            } while(true);
            is.clear();
            poscar.numOfEachType.squeeze();
        }
        /* Read format type */ {
            const int ch = std::tolower(is.get());
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (ch == 'd')
                poscar.type = Poscar::Direct;
            else if (ch == 'c')
                poscar.type = Poscar::Cartesian;
            else
                throw BadFileFormatException();
        }
        /* Read atom pos */ {
            const size_t atomCount = poscar.getAtomCount();
            poscar.pos.resize(atomCount, 3);
            size_t i = 0;
            for (; i < atomCount - 1; i++) {
                is >> poscar.pos(i, 0) >> poscar.pos(i, 1) >> poscar.pos(i, 2);
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            is >> poscar.pos(i, 0) >> poscar.pos(i, 1) >> poscar.pos(i, 2);
        }
        if (!is)
            throw BadFileFormatException();
        return is;
    }

    void Poscar::standrizeLattice() {
        using MatrixType = typename LatticeMatrix::ColMatrix;
        MatrixType temp = lattice.transpose();
        using VectorType = Vector<ScalarType, 3>;
        VectorType buffer{};
        householder(temp.col(0), buffer);
        applyHouseholder(buffer, temp);

        auto buffer1 = buffer.head(2);
        householder(temp.col(1).tail(1), buffer1);
        auto corner = temp.bottomRightCorner(1);
        applyHouseholder(buffer1, corner);
        lattice = temp.transpose();
    }
    /**
     * Extend the cell in z direction, with all distance of atoms in cell not changed.
     */
    void Poscar::extendInZ(ScalarType factor) {
        assert(type == Direct);
        assert(lattice(0, 1).isZero());
        assert(lattice(0, 2).isZero());
        assert(lattice(1, 2).isZero());
        lattice(2, 2) *= factor;

        const ScalarType inv_factor = reciprocal(factor);
        auto col = pos.col(2);
        col.asVector() *= inv_factor;

        const ScalarType temp = ScalarType::One() - inv_factor;
        for (size_t i = 0; i < col.getLength(); ++i)
            if (col[i] > ScalarType(0.5))
                col[i] += temp;
    }

    size_t Poscar::getAtomCount() const noexcept {
        size_t result = 0;
        for (size_t i = 0; i < numOfEachType.getLength(); ++i)
            result += numOfEachType[i];
        return result;
    }
}
