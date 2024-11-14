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

#include <fstream>
#include <Physica/Core/MultiPrecision/Real.h>

namespace Physica::Core {
    class PHYSICA_API Potcar {
        using ScalarType = float32;
        ScalarType mass;
        unsigned char numValenceElectron;
    public:
        explicit Potcar(const char* path);
        Potcar(const Potcar&) = default;
        Potcar(Potcar&&) noexcept = default;
        ~Potcar() = default;
        /* Operators */
        Potcar& operator=(Potcar obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(Potcar& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getMass() const noexcept { return mass; }
        [[nodiscard]] unsigned int getNumValenceElectron() const noexcept { return numValenceElectron; }
    private:
        void readFile(std::ifstream& fin);
    };
}
