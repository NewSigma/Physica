/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Parallel/PBSWrapper.h"
#include <cerrno>
#include <fstream>
#include "Physica/Core/Exception/BadFileFormatException.h"

using namespace Physica;

const PBSWrapper& PBSWrapper::getInstance() noexcept {
    static PBSWrapper pbs{};
    return pbs;
}

PBSWrapper::PBSWrapper() : jobCore(0) {
    readJobCore();
    readHostList();
}

void PBSWrapper::readJobCore() {
    char* jobCoreStr = getenv("PBS_NP");
    if (jobCoreStr == nullptr)
        return;

    errno = 0;
    jobCore = strtoul(jobCoreStr, nullptr, 10);
    const bool isOverFlow = errno != 0;
    if (isOverFlow)
        jobCore = 0;
}

void PBSWrapper::readHostList() {
    char* path_nodefile = getenv("PBS_NODEFILE");
    if (path_nodefile == nullptr)
        return;
    hostList.resize(jobCore);

    char buffer[hostLength];
    std::ifstream fin(path_nodefile);
    for (unsigned int i = 0; i < jobCore; ++i) {
        fin.getline(buffer, hostLength);
        hostList[i] = std::string(buffer);
    }

    if (!fin)
        throw BadFileFormatException("PBS_NODEFILE");
}
