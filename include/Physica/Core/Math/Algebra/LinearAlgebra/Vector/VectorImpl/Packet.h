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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
    #include "vectorclass/vectorclass.h"
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

namespace Physica::Core::Internal {
    template<>
    class Traits<Vec4f> {
    public:
        using ScalarType = float;
    };

    template<>
    class Traits<Vec8f> {
    public:
        using ScalarType = float;
    };

    template<>
    class Traits<Vec16f> {
    public:
        using ScalarType = float;
    };

    template<>
    class Traits<Vec2d> {
    public:
        using ScalarType = double;
    };

    template<>
    class Traits<Vec4d> {
    public:
        using ScalarType = double;
    };

    template<>
    class Traits<Vec8d> {
    public:
        using ScalarType = double;
    };
}
