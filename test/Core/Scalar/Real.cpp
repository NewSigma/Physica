/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Scalar/Real.h"

using namespace Physica;

int main() {
    // We handle type casts in the CRTP base class. Test that we do not accidentally fall into infinite recursion.
    std::ignore = float32(0) < float64(0);
    std::ignore = float32(0) > float64(0);
    std::ignore = float32(0) <= float64(0);
    std::ignore = float32(0) >= float64(0);
    return 0;
}
