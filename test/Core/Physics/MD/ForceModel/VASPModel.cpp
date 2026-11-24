/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
