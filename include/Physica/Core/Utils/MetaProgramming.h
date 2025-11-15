/*
 * Copyright 2024-2025 Weibo He.
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

namespace Physica {
    template<class> class ScalarRef;

    namespace Internal {
        template<class T>
        struct remove_ScalarRef {
            using Type = T;
        };

        template<class T>
        struct remove_ScalarRef<ScalarRef<T>> {
            using Type = T;
        };

        template<template<class...> class, class>
        struct is_instance_of : std::false_type {};

        template<template<class...> class Template, class... Args>
        struct is_instance_of<Template, Template<Args...>> : std::true_type {};

        template<template<auto, class...> class, class>
        struct is_instance_of_xt : std::false_type {};

        template<template<auto, class...> class Template, auto Arg0, class... Args>
        struct is_instance_of_xt<Template, Template<Arg0, Args...>> : std::true_type {};

        template<template<class, auto...> class, class>
        struct is_instance_of_tx : std::false_type {};

        template<template<class, auto...> class Template, class Arg0, auto... Args>
        struct is_instance_of_tx<Template, Template<Arg0, Args...>> : std::true_type {};
    }

    template<class T>
    struct remove_cvref {
        using Type = Internal::remove_ScalarRef<std::remove_cvref_t<T>>::Type;
    };

    template<template<class...> class Template, class T>
    concept instanceof = Internal::is_instance_of<Template, std::remove_cvref_t<T>>::value;
    /**
     * x: Non type template param
     * t: Type template param
     *
     * It is unfortunate that we cannot implement instanceof in an elegant way.
     */
    template<template<auto, class...> class Template, class T>
    concept instanceof_xt = Internal::is_instance_of_xt<Template, std::remove_cvref_t<T>>::value;

    template<template<class, auto...> class Template, class T>
    concept instanceof_tx = Internal::is_instance_of_tx<Template, std::remove_cvref_t<T>>::value;
}
