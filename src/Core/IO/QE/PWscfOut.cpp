/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/IO/QE/PWscfOut.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Physics/PhyConst.h"

using namespace Physica;

PWscfOut::PWscfOut(const char* path, size_t numAtom)
        : fin(path), force(3 * numAtom), buffer(DefaultBufferSize) {
    if (!fin)
        throw IOException("[Error]: File not found");
    readForce();
}

VectorND<typename PWscfOut::ScalarType> PWscfOut::makeTotalForces() {
    fin.seekg(std::ios::beg);

    VectorND<ScalarType> result{};
    std::string str{};
    do {
        fin.getline(buffer.data(), buffer.getLength());
        if (!fin) [[unlikely]]
            break;
        str = buffer.data();
        const bool success = str.find("Total force") != std::string::npos;
        if (success) {
            /* Back */
            const size_t count = fin.gcount();
            for (size_t i = 0; i < count; ++i)
                fin.unget();
            /* Forward */
            fin.getline(buffer.data(), buffer.getLength(), '=');
            ScalarType temp{};
            fin >> temp;
            result.append(temp);
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    } while (fin.good());

    if (!fin.eof())
        throw BadFileFormatException("[Error]: Failed to read force");

    constexpr double convertFactor = PhyConst<AU>::planck / PhyConst<QE>::planck;
    force *= ScalarType(convertFactor);
    return result;
}

void PWscfOut::swap(PWscfOut& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    fin.swap(obj.fin);
    force.swap(obj.force);
}

void PWscfOut::readForce() {
    std::string str{};
    do {
        fin.getline(buffer.data(), buffer.getLength());
        str = buffer.data();
        assert(str.size() < DefaultBufferSize && "[Error]: Unexpected log length per line");
        const bool success = str.find("Forces acting on atoms") != std::string::npos;
        if (success) {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (size_t atom = 0; atom < getNumAtom(); ++atom) {
                fin.getline(buffer.data(), buffer.getLength(), '=');
                for (int dim = 0; dim < 3; ++dim)
                    fin >> force[atom * 3 + dim];
            }
            break;
        }
    } while (fin.good());

    if (!fin)
        throw BadFileFormatException("[Error]: Failed to read force");

    constexpr double convertFactor = PhyConst<AU>::planck / PhyConst<QE>::planck;
    force *= ScalarType(convertFactor);
}
