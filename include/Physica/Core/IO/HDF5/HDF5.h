/*
 * Copyright 2023 Weibo He.
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
#pragma once

#include "Physica/Macro.h"

#ifdef PHYSICA_HDF5
    #include <H5Cpp.h>

    #include "H5DataSet.h"
    #include "H5DataSpace.h"
    #include "H5Location.h"
    #include "H5File.h"
    #include "H5Group.h"
#else
    namespace H5 {
        class DataType;
        class DSetMemXferPropList;
    }

    namespace Physica::Core {
        template<class Derived> class DataSpaceBase {};
        template<size_t Dim> class H5DataSpace {};
        template<size_t Dim> class H5DataSet {};
        class H5File;
        class H5Group {};
        class H5Location {};
    }
#endif
