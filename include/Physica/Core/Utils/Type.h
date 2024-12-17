/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<class ScalarType>
    class ScalarRef;

    namespace Internal {
        template<class T>
        struct remove_ScalarRef {
            using Type = T;
        };

        template<class T>
        struct remove_ScalarRef<ScalarRef<T>> {
            using Type = T;
        };
    }

    template<class T>
    struct remove_cvref {
        using Type = Internal::remove_ScalarRef<std::remove_cvref_t<T>>::Type;
    };
}
