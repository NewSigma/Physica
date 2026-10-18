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
#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Physica/Core/IO/HDF5/H5File.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using PhononType = FrozenPhonon<T>;
using GridType = PhononType::RSpaceFCGrid;
using RandomSource = Random<>;

namespace {
    void hdf5() {
        TempFile temp("/tmp/tmpXXXXXX");
        const size_t dof = 6;
        const auto grid = GridType::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ, dof, dof);
        {
            auto h5f = H5File::open(temp.getName());
            PhononType::write(grid, h5f, "Si");
        }
        {
            auto h5f = H5File::open(temp.getName(), H5File::ReadOnly);
            GridType read{};
            PhononType::read(read, h5f, "Si");
            expect(read.getShape() == grid.getShape());
            expect(read == grid);
        }
    }
}

int main() {
    hdf5();
    return 0;
}
