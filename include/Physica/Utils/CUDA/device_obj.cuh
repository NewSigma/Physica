/*
 * Copyright 2022-2023 WeiBo He.
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
 * \class device_obj: Provide non-invasive implementation for device, which is determined by the nature of nvcc seperate compiling.
 * 
 * Resource is maintained by host and can be used on device.
 * 
 * Class name is compatible to \class thrust::device_ptr and \class thrust::device_reference.
 */
namespace Physica {
    namespace Core {
        namespace Internal {
            template<class T> class device_obj;
        }

        template<class T> class device_obj;
    }

    namespace Utils {
        namespace Internal {
            template<class T, bool IsTrivial>
            struct add_device_obj_impl {
                using Type = T;
            };

            template<class T>
            struct add_device_obj_impl<T, false> {
                using Type = typename T::device_obj_type;
            };
        }

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
        struct is_device_obj<Physica::Core::device_obj<T>> {
            constexpr static bool value = true;
        };

        template<class T>
        struct is_device_obj<Physica::Core::Internal::device_obj<T>> {
            constexpr static bool value = true;
        };

        template<class T>
        struct add_device_obj {
            using Type = typename Internal::add_device_obj_impl<T, std::is_trivial<T>::value>::Type;
        };
    }
}
