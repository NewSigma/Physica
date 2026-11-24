/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/SolidState/VASP/Potcar.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Exception/IOException.h"

using namespace Physica;

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
    constexpr size_t DefaultBufferSize = 1024; // 1024 shall be enough
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
    } while (fin.good());

    if (!fin)
        throw BadFileFormatException("[Error]: Failed to read POTCAR");
}
