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
#include <fstream>
#include "Physica/Core/Physics/MD/ForceModel/QEModel.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using QEModelType = QEModel<T>;
using ElementTypeArray = Poscar<T>::ElementTypeArray;

namespace {
    void constructor() {
        TempFile input("/tmp/tmpXXXXXX");
        ElementTypeArray elementTypes{8, 8, 16};
        const QEModelType model("pw.x", input.getName(), elementTypes, 4);

        expect(model.getNumParticle() == 3);
        expect(model.getNumMPIProcess() == 4);
    }
}

int main() {
    constructor();
    return 0;
}
