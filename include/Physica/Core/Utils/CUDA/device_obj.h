/*
 * Copyright 2022-2026 Weibo He.
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

    namespace Internal {
        template<class T>
        struct remove_device_obj_trivial {
            static_assert(!std::is_reference_v<T> && !std::is_const_v<T> && !is_codiff_v<T>, "[Error]: Not a trivial type");
            using type = T;
        };

        template<class T>
        struct remove_device_obj_trivial<device_obj<T>> {
            using type = T;
        };
    }

    template<class T>
    struct is_device_obj {
        constexpr static bool value = instanceof<device_obj, remove_codiff_t<std::remove_cvref_t<T>>>;
    };

    template<class T>
    constexpr bool is_device_obj_v = is_device_obj<T>::value;

    template<class T>
    struct remove_device_obj {
    private:
        template<bool NoCVRef>
        struct Helper {
            using type = Internal::remove_device_obj_trivial<T>::type;
        };

        template<>
        struct Helper<false> {
            using U1 = Internal::remove_device_obj_trivial<remove_codiff_t<std::remove_cvref_t<T>>>::type;
            using U2 = std::conditional<is_codiff<T>::value, CoDiff<U1>, U1>::type;
            using type = copy_cvref<T, U2>::type;
        };
    public:
        using type = Helper<!std::is_reference_v<T> && !std::is_const_v<T>>::type;
    };

    template<class T>
    using remove_device_obj_t = remove_device_obj<T>::type;

    template<class T>
    struct add_device_obj {
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj is not allowed");
    private:
        template<bool Trivial>
        struct Helper {
            using type = device_obj<T>;
        };

        template<>
        struct Helper<false> {
            using U1 = device_obj<remove_codiff_t<std::remove_cvref_t<T>>>;
            using U2 = std::conditional<is_codiff<T>::value, CoDiff<U1>, U1>::type;
            using type = copy_cvref<T, U2>::type;
        };
    public:
        using type = Helper<!std::is_reference_v<T> && !std::is_const_v<T> && !is_codiff_v<T>>::type;
    };

    template<class T>
    using add_device_obj_t = add_device_obj<T>::type;

    template<class T>
    concept CUDA = is_device_obj<T>::value;

    template<CUDA T>
    class device_obj<T> {
        static_assert(!CUDA<T>, "[Error]: Nested device_obj is not allowed");
    };
}
