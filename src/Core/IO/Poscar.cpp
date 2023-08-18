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
#include <unordered_set>
#include <algorithm>
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Exception/NotImplementedException.h"
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    Poscar::Poscar()
            : Base(), elementTypes(), numOfEachType() {}

    Poscar::Poscar(Base base, ElementTypeArray elementTypes_, Utils::Array<size_t> numOfEachType_)
            : Base(std::move(base))
            , elementTypes(std::move(elementTypes_))
            , numOfEachType(std::move(numOfEachType_)) {
        assert(getAtomCount() == sumNumOfEachType());
        assert(elementTypes.getLength() == numOfEachType.getLength());
    }

    Poscar::Poscar(CrystalCell cell) : Base(std::move(cell)) {
        /* Get size */ {
            std::unordered_set<uint16_t> set{};
            for (uint16_t elem : cell.getAtomicNumbers())
                set.insert(elem);
            const size_t size = set.size();
            elementTypes.resize(size);
            numOfEachType.resize(size, 0);
        }
        size_t index = 0;
        uint16_t temp = cell.getAtomicNumbers()[0];
        for (auto atomicNumber : cell.getAtomicNumbers()) {
            const bool isSame = temp == atomicNumber;
            if (!isSame) {
                index += 1;
                elementTypes[index] = atomicNumber;
            }
            numOfEachType[index] += 1;
            temp = atomicNumber;
        }
    }

    std::ostream& operator<<(std::ostream& os, const Poscar& poscar) {
        os << '\n';
        os << 1.0 << '\n';
        os << poscar.lattice;
        if (!poscar.elementTypes.empty()) {
            for (auto type : poscar.elementTypes)
                os << ' ' << PhyConst<SI>::elementSymbol[type];
            os << '\n';
        }
        for (size_t i = 0; i < poscar.numOfEachType.getLength(); ++i)
            os << ' ' << poscar.numOfEachType[i];
        os << '\n';
        os << ((poscar.type == Poscar::Type::Direct) ? "Direct\n" : "Cartesian\n");
        os << poscar.pos;
        return os;
    }

    std::istream& operator>>(std::istream& is, Poscar& poscar) {
        assert(is.good());
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        is >> poscar.lattice(0, 0) >> poscar.lattice(0, 1) >> poscar.lattice(0, 2);
        is >> poscar.lattice(1, 0) >> poscar.lattice(1, 1) >> poscar.lattice(1, 2);
        is >> poscar.lattice(2, 0) >> poscar.lattice(2, 1) >> poscar.lattice(2, 2);

        poscar.readTypesAndNumber(is);
        /* Read format type */ {
            const int ch = std::tolower(is.get());
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (ch == 'd')
                poscar.type = Poscar::Type::Direct;
            else if (ch == 'c')
                poscar.type = Poscar::Type::Cartesian;
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
        temp(1, 0) = temp(2, 0) = temp(2, 1) = ScalarType(0); //Clear numeric error
        lattice = temp.transpose();
    }
    /**
     * Extend the cell in z direction, with all distance of atoms in cell not changed. Use for 2D material only
     */
    void Poscar::extendInZ(ScalarType factor) {
        if (Base::type == Type::Cartesian) {
            Base::toDirect();
            extendInZ_direct(factor);
            Base::toCartesian();
        }
        else
            extendInZ_direct(factor);
    }

    void Poscar::toUnitCell(unsigned int x, unsigned int y, unsigned int z) {
        Base::toUnitCell(x, y, z);
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

    void Poscar::swap(Poscar& poscar) noexcept {
        Base::swap(poscar);
        elementTypes.swap(poscar.elementTypes);
        numOfEachType.swap(poscar.numOfEachType);
    }

    void Poscar::readTypesAndNumber(std::istream& is) {
        size_t count = 0;
        /* Read types */ {
            is >> count;
            const bool hasSymbol = is.fail();
            if (hasSymbol) {
                is.clear();

                elementTypes.reserve(8); // 8 is enough for common use
                int ch = is.peek();
                while (true) {
                    while (std::isspace(ch) && ch != '\n') {
                        is.get();
                        ch = is.peek();
                    }
                    if (ch == '\n')
                        break;
                    is.get();

                    const int next_ch = is.peek();
                    const bool isAlpha = std::isalpha(next_ch);
                    const int ch1 = isAlpha ? next_ch : '\0';
                    elementTypes.append(elementSymbolToNumber(ch, ch1));
                    if (isAlpha) {
                        is.get();
                        ch = is.peek();
                    }
                    else
                        ch = next_ch;
                }
                elementTypes.squeeze();
                is.get();
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
        const size_t atomCount = sumNumOfEachType();
        pos.resize(atomCount, 3);
        size_t i = 0;
        for (; i < atomCount - 1; i++) {
            is >> pos(i, 0) >> pos(i, 1) >> pos(i, 2);
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        is >> pos(i, 0) >> pos(i, 1) >> pos(i, 2);
    }

    size_t Poscar::sumNumOfEachType() const {
        size_t result = 0;
        for (size_t num : numOfEachType)
            result += num;
        return result;
    }

    void Poscar::extendInZ_direct(ScalarType factor) {
        assert(type == Type::Direct);
        assert(lattice(0, 2).isZero());
        assert(lattice(1, 2).isZero());
        lattice(2, 2) *= factor;

        const ScalarType inv_factor = Core::reciprocal(factor);

        auto col = pos.col(2);
        for (size_t i = 0; i < col.getLength(); ++i) {
            if (col[i] > ScalarType(0.5))
                col[i] += (ScalarType(1) - col[i]) * inv_factor;
            else
                col[i] *= inv_factor;
        }
    }

    uint8_t Poscar::elementSymbolToNumber(char ch1, char ch2) {
        assert(std::isalpha(ch1) && (std::isalpha(ch2) || ch2 == '\0'));
        uint8_t i = 0;
        for (const char* symbol : PhyConst<SI>::elementSymbol) {
            if (symbol[0] == ch1 && symbol[1] == ch2)
                return i;
            i += 1;
        }
        throw NotImplementedException("Unrecognized element symbol");
    }
}
