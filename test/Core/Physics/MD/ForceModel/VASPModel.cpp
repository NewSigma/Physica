/*
 * Copyright 2026 Weibo He.
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
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "Physica/Core/Parallel/Algorithm/Sequential.h"
#include "Physica/Core/Physics/MD/ForceModel/VASPModel.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using VASPModelType = VASPModel<T>;

namespace {
    void constructor() {
        TempFile incar("/tmp/tmpXXXXXX");
        TempFile potcar("/tmp/tmpXXXXXX");
        TempFile kpoints("/tmp/tmpXXXXXX");
        const VASPModelType model("vasp", incar.getName(), potcar.getName(), kpoints.getName(), Array<size_t>{1, 2}, 4);

        expect(model.getNumParticle() == 3);
        const std::filesystem::path workingDir = model.getWorkingDir().getName();
        expect(std::filesystem::exists(workingDir / "log"));
    }
}

int main() {
    constructor();
    return 0;
}
