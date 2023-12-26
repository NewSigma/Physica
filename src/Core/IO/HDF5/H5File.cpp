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
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Exception/IOException.h"

namespace Physica::Core {
    H5File::H5File(
            const char* name,
            unsigned int openflag_,
            const H5::FileCreatPropList& create_plist,
            const H5::FileAccPropList& access_plist) : Base(name, openflag_, create_plist, access_plist), openflag(openflag_) {}
    
    H5File::H5File(const H5File& obj) : Base(obj), openflag(obj.openflag) {}

    H5File& H5File::operator=(H5File& obj) {
        Base::operator=(obj);
        openflag = obj.openflag;
        return *this;
    }

    H5DataSet<1> H5File::createDataSet(const char* filepath, const char* name) {
        std::ifstream fin(filepath);
        if (!fin)
            throw IOException("[Error]: File not found");
        fin.seekg(0, std::ios::end);
        const hsize_t size = fin.tellg();
        fin.seekg(0, std::ios::beg);

        auto buffer = Utils::Array<char>(size);
        fin.read(buffer.data(), size);
        const auto space = H5DataSpace<1>({size});
        auto dataset = Base::createDataSet(name, H5::PredType::NATIVE_CHAR, space);
        dataset.write(buffer.data(), H5::PredType::NATIVE_CHAR);
        return dataset;
    }
}
