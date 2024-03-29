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
#include "Physica/Utils/TestHelper.h"

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

        poscar.readNumOfEachType(is);
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
        poscar.readAtomPos(is);
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
        for (int i = 0; i < 3; ++i) {
            if (temp(i, i).isNegative()) {
                auto row = temp.row(i);
                row = -row.asVector();
            }
        }
        temp(1, 0) = temp(2, 0) = temp(2, 1) = ScalarType::Zero(); //Clear numeric error
        lattice = temp.transpose();
    }
    /**
     * Extend the cell in z direction, with all distance of atoms in cell not changed. Use for 2D material only
     */
    void Poscar::extendInZ(ScalarType factor) {
        assert(type == Direct);
        assert(lattice(0, 2).isZero());
        assert(lattice(1, 2).isZero());
        lattice(2, 2) *= factor;

        const ScalarType inv_factor = reciprocal(factor);

        auto col = pos.col(2);
        for (size_t i = 0; i < col.getLength(); ++i) {
            if (col[i] > ScalarType(0.5))
                col[i] += (ScalarType::One() - col[i]) * inv_factor;
            else
                col[i] *= inv_factor;
        }
    }
    /**
     * Helper for phonopy, convert from super cell to unit cell
     */
    void Poscar::superToUnit(size_t x, size_t y, size_t z) {
        assert(x > 0 && y > 0 && z > 0);
        assert(pos.getRow() % (x * y * z) == 0);
        const ScalarType inv_x = reciprocal(ScalarType(x));
        const ScalarType inv_y = reciprocal(ScalarType(y));
        const ScalarType inv_z = reciprocal(ScalarType(z));

        auto rowX = lattice.row(0);
        rowX *= inv_x;
        auto rowY = lattice.row(1);
        rowY *= inv_y;
        auto rowZ = lattice.row(2);
        rowZ *= inv_z;

        auto colX = pos.col(0);
        colX *= ScalarType(x);
        auto colY = pos.col(1);
        colY *= ScalarType(y);
        auto colZ = pos.col(2);
        colZ *= ScalarType(z);

        const size_t newPosSize = pos.getRow() / (x * y * z);
        PositionMatrix newPos(newPosSize, 3);
        size_t toFill = 0;
        size_t toCheck = 0;
        const ScalarType one = ScalarType::One();
        for (; toFill < newPosSize; ++toFill) {
            for (; toCheck < pos.getRow(); ++toCheck) {
                auto rowToCheck = pos.row(toCheck);
                if (rowToCheck[0] <= one && rowToCheck[1] <= one && rowToCheck[2] <= one) {
                    auto rowToFill = newPos.row(toFill);
                    rowToFill = rowToCheck.asVector();
                    ++toCheck;
                    break;
                }
            }
        }
        pos = newPos;

        for (size_t& num : numOfEachType)
            num /= (x * y * z);
    }

    Poscar::CrystalSystem Poscar::getCrystalSystem(double precision) const noexcept {
        using namespace Utils;
        const ScalarType norm_list[3]{lattice.row(0).squaredNorm(),
                                      lattice.row(1).squaredNorm(),
                                      lattice.row(2).squaredNorm()};
        const ScalarType angle_list[3]{lattice.row(1).angleTo(lattice.row(2)),
                                       lattice.row(0).angleTo(lattice.row(2)),
                                       lattice.row(0).angleTo(lattice.row(1))};
        const bool norm_equal_list[3]{scalarNear(norm_list[0], norm_list[1], precision),
                                      scalarNear(norm_list[1], norm_list[2], precision),
                                      scalarNear(norm_list[2], norm_list[0], precision)};
        const bool angle_equal_list[3]{scalarNear(angle_list[0], angle_list[1], precision),
                                       scalarNear(angle_list[1], angle_list[2], precision),
                                       scalarNear(angle_list[2], angle_list[0], precision)};

        const bool allAngleSame = angle_equal_list[0] && angle_equal_list[1];
        if (allAngleSame) {
            if (scalarNear(angle_list[0], ScalarType(M_PI_2), precision)) {
                const unsigned int sameNormCount = norm_equal_list[0] + norm_equal_list[1] + norm_equal_list[2];
                [[maybe_unused]] const bool compareTransitive = sameNormCount != 2;
                assert(compareTransitive);
                if (sameNormCount == 3)
                    return Cubic;
                if (sameNormCount == 1)
                    return Tetragonal;
                return Orthohombic;
            }
            return Rhombohedral;
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (angle_equal_list[i]) {
                    if (scalarNear(angle_list[i], ScalarType(M_PI_2), precision)) {
                        if (scalarNear(angle_list[(i + 2) % 3], ScalarType(M_PI * 2 / 3), precision))
                            return Hexagonal;
                        return Monoclinic;
                    }
                    [[maybe_unused]] const bool onlyOneSamePair = !angle_equal_list[(i + 1) % 3] && !angle_equal_list[(i + 2) % 3];
                    assert(onlyOneSamePair);
                    break;
                }
            }
        }
        return Triclinic;
    }

    size_t Poscar::getAtomCount() const noexcept {
        size_t result = 0;
        for (size_t i = 0; i < numOfEachType.getLength(); ++i)
            result += numOfEachType[i];
        return result;
    }

    void Poscar::swap(Poscar& poscar) noexcept {
        lattice.swap(poscar.lattice);
        pos.swap(poscar.pos);
        numOfEachType.swap(poscar.numOfEachType);
        std::swap(type, poscar.type);
    }

    void Poscar::readNumOfEachType(std::istream& is) {
        size_t count = 0;
        /* Get first type count */ {
            is >> count;
            if (is.fail()) {
                is.clear();
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                is >> count;
            }

            if (!is)
                throw BadFileFormatException();
        }
        numOfEachType.reserve(8);
        do {
            numOfEachType << count;
            is >> count;
        } while(is.good());
        is.clear();
        numOfEachType.squeeze();
    }

    void Poscar::readAtomPos(std::istream& is) {
        const size_t atomCount = getAtomCount();
        pos.resize(atomCount, 3);
        size_t i = 0;
        for (; i < atomCount - 1; i++) {
            is >> pos(i, 0) >> pos(i, 1) >> pos(i, 2);
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        is >> pos(i, 0) >> pos(i, 1) >> pos(i, 2);
    }
}
