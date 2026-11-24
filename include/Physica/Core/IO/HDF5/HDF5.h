/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"

#ifdef PHYSICA_HDF5
    #include <hdf5.h>

    #include "H5Type.h"
    #include "H5Dataset.h"
    #include "H5DataSpace.h"
    #include "H5File.h"
    #include "H5Group.h"
#else
    namespace Physica {
        template<class Derived> class DataSpaceMixin {};
        template<size_t Dim> class H5DataSpace {};
        template<size_t Dim> class H5Dataset {};
        class H5Attribute;
        class H5DataType;
        class H5File;
        class H5Group {};
        class H5Loc {};
    }
#endif
