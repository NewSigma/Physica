/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/IO/VASP/Potcar.h"
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica {
    Potcar::Potcar(const char* path) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: No POTCAR file found");
        readFile(fin);
    }

    void Potcar::swap(Potcar& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mass.swap(obj.mass);
        std::swap(numValenceElectron, obj.numValenceElectron);
    }

    void Potcar::readFile(std::ifstream& fin) {
        constexpr size_t DefaultBufferSize = 1024; //1024 shall be enough
        Array<char> buffer(DefaultBufferSize);
        std::string str{};
        do {
            fin.getline(buffer.data(), buffer.getLength());
            if (!fin) [[unlikely]]
                break;
            str = buffer.data();
            const bool success = str.find("ZVAL") != std::string::npos;
            if (success) {
                float pomass;
                sscanf(buffer.data(), " POMASS = %f; ZVAL = %hhd", &pomass, &numValenceElectron);
                mass = ScalarType(pomass);
                break;
            }
        } while(fin.good());

        if (!fin)
            throw BadFileFormatException("[Error]: Failed to read POTCAR");
    }
}
