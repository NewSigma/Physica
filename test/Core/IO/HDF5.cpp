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
#include <iostream>
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Utils;

int main() {
    TempFile temp("/tmp/tmpXXXXXX");
    const char* str = "This is a str\nstr line2";
    
    const auto dataspace = H5DataSpace<1>({strlen(str)});
    {
        H5File h5f(temp.getName(), H5File::ReadWrite | H5File::Creat);
        auto dataset = h5f.createDataSet("/data", H5::PredType::NATIVE_CHAR, dataspace);
        dataset.write(str, H5::PredType::NATIVE_CHAR);
    }
    {
        char buffer[32];
        H5File h5f(temp.getName(), H5File::ReadOnly);
        auto dataset = h5f.openDataSet<1>("/data");
        dataset.readStr(buffer);
        if (strcmp(str, buffer) != 0)
            return 1;
    }
    return 0;
}
