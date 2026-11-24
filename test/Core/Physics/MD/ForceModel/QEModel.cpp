/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
