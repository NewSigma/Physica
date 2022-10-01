/*
 * Copyright 2022 WeiBo He.
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
/**
 * \class device_obj: Provide non-invasive implementation for device
 * 
 * Class name is compatible to \class thrust::device_ptr and \class thrust::device_reference.
 */
namespace Physica::Core {
    namespace Internal {
        template<class T> class device_obj final : public T {
            static_assert(std::is_trivial<T>::value, "[Error]: non-trivial class does not have default device_obj");
        };
    }
    template<class T> class device_obj final : public T {
        static_assert(std::is_trivial<T>::value, "[Error]: non-trivial class does not have default device_obj");
    };
}

namespace Physica::Utils {

    template<class T> class device_obj final : public T {
        static_assert(std::is_trivial<T>::value, "[Error]: non-trivial class does not have default device_obj");
    };

    template<class T>
    struct remove_device_obj {
        using Type = T;
    };

    template<class T>
    struct remove_device_obj<device_obj<T>> {
        using Type = T;
    };
}
