/*
 * Copyright 2022-2025 Weibo He.
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

#include <type_traits>
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    /**
     * \class device_obj: Provide non-invasive implementation for device, which is determined by the nature of CUDA seperate compiling.
     * 
     * Resource is maintained by host and can be used on device.
     * 
     * Class name is compatible to \class thrust::device_ptr and \class thrust::device_reference.
     */
    template<class T> class device_obj;

    template<class T>
    struct is_device_obj {
        constexpr static bool value = false;
    };

    template<class T>
    struct is_device_obj<device_obj<T>> {
        constexpr static bool value = true;
    };

    template<class T>
    struct remove_device_obj {
        using Type = T;
    };

    template<class T>
    struct remove_device_obj<device_obj<T>> {
        using Type = T;
    };

    template<class T>
    using remove_device_obj_t = remove_device_obj<T>::Type;

    template<class T>
    concept CUDA = is_device_obj<remove_codiff_t<std::remove_cvref_t<T>>>::value;

    template<CUDA T>
    class device_obj<T> {
        static_assert(!CUDA<T>, "[Error]: Nested device_obj is not allowed");
    };
}
