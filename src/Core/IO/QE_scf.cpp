/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/IO/QE_scf.h"
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    QE_scf::QE_scf(const char* path, size_t numAtom) : force(3 * numAtom) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: No QE output file found");
        fin.seekg(0, std::ios::end);
        const auto size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        Utils::Array<char> buffer(size);

        readForce(fin, buffer);
    }

    QE_scf& QE_scf::operator=(QE_scf obj) noexcept {
        swap(obj);
        return *this;
    }

    void QE_scf::swap(QE_scf& obj) noexcept {
        force.swap(obj.force);
    }

    void QE_scf::readForce(std::ifstream& fin, Utils::Array<char>& buffer) {
        std::string str{};
        do {
            fin.getline(buffer.data(), buffer.getLength());
            str = buffer.data();
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
        } while(bool(fin));

        if (!fin)
            throw BadFileFormatException("[Error]: Failed to read force");

        constexpr double convertFactor = PhyConst<AU>::planck / PhyConst<QE>::planck;
        force *= ScalarType(convertFactor);
    }
}
