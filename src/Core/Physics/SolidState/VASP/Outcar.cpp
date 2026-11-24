/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/SolidState/VASP/Outcar.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/PhyConst.h"

using namespace Physica;

Outcar::Outcar(const char* path, unsigned int numAtom)
        : force(numAtom * 3) {
    std::ifstream fin(path);
    if (!fin)
        throw IOException("[Error]: No OUTCAR file found");
    fin.seekg(0, std::ios::end);
    const auto size = fin.tellg();
    fin.seekg(0, std::ios::beg);
    Array<char> buffer(size);

    readForce(fin, buffer);
    readEnergy(fin, buffer);
}

void Outcar::swap(This& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    force.swap(obj.force);
    internalEnergy.swap(obj.internalEnergy);
}

void Outcar::readForce(std::ifstream& fin, Array<char>& buffer) {
    using MatrixType = DenseMatrix<ScalarType, MatrixMajor::Row, Dynamic, 6>;

    std::string str{};
    do {
        fin.getline(buffer.data(), buffer.getLength());
        str = buffer.data();
        const bool success = str.find("TOTAL-FORCE") != std::string::npos;
        if (success) {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            MatrixType pos_force(getNumAtom(), 6);
            for (auto& elem : pos_force.asArray())
                fin >> elem;

            size_t index = 0;
            for (size_t r = 0; r < pos_force.getRow(); ++r) {
                for (size_t c = 3; c < pos_force.getCol(); ++c) {
                    force[index] = pos_force[r, c];
                    ++index;
                }
            }
            break;
        }
    } while (fin.good());

    if (!fin)
        throw BadFileFormatException("[Error]: Failed to read force");
    force *= ScalarType(PhyConst<AU>::eVToHartree(1) * PhyConst<AU>::bohrToAngstrom(1));
}

void Outcar::readEnergy(std::ifstream& fin, Array<char>& buffer) {
    std::string str{};
    do {
        fin.getline(buffer.data(), buffer.getLength());
        str = buffer.data();
        const bool success = str.find("energy(sigma->0)") != std::string::npos;
        if (success) {
            unsigned int count = 0;
            for (; count < str.size(); ++count)
                if (str[count] == '=')
                    break;
            ++count;
            str.erase(0, count);
            internalEnergy = std::stod(str, nullptr);
            break;
        }
    } while (fin.good());

    if (!fin)
        throw BadFileFormatException("[Error]: Failed to read energy");
    internalEnergy *= ScalarType(PhyConst<AU>::eVToHartree(1));
}
