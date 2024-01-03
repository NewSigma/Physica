/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/IO/Outcar.h"
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    Outcar::Outcar(const char* path, unsigned int numAtom) : force(numAtom * 3) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: No OUTCAR file found");
        fin.seekg(0, std::ios::end);
        const auto size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        Utils::Array<char> buffer(size);

        readForce(fin, buffer);
        readEnergy(fin, buffer);
    }

    void Outcar::swap(Outcar& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        force.swap(obj.force);
        internalEnergy.swap(obj.internalEnergy);
    }

    void Outcar::readForce(std::ifstream& fin, Utils::Array<char>& buffer) {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dynamic, 6>;

        std::string str{};
        do {
            fin.getline(buffer.data(), buffer.getLength());
            str = buffer.data();
            const bool success = str.find("TOTAL-FORCE") != std::string::npos;
            if (success) {
                fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                MatrixType pos_force(getNumAtom(), 6);
                fin >> pos_force;
                size_t index = 0;
                for (size_t r = 0; r < pos_force.getRow(); ++r) {
                    for (size_t c = 3; c < pos_force.getColumn(); ++c) {
                        force[index] = pos_force(r, c);
                        ++index;
                    }
                }
                break;
            }
        } while(fin.good());

        if (!fin)
            throw BadFileFormatException("[Error]: Failed to read force");
    }

    void Outcar::readEnergy(std::ifstream& fin, Utils::Array<char>& buffer) {
        std::string str{};
        do {
            fin.getline(buffer.data(), buffer.getLength());
            str = buffer.data();
            const bool success = str.find("energy(sigma->0)") != std::string::npos;
            if (success) {
                unsigned int count = 0;
                for (;count < str.size(); ++count)
                    if (str[count] == '=')
                        break;
                ++count;
                str.erase(0, count);
                internalEnergy = std::stod(str, nullptr);
                break;
            }
        } while(fin.good());

        if (!fin)
            throw BadFileFormatException("[Error]: Failed to read energy");
    }
}
