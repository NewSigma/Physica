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

#include <fstream>
#include <utility>
#include "Poscar.h"

namespace Physica::Core {
    class PHYSICA_API Xdatcar {
        using ScalarType = Scalar<Double>;
    private:
        Poscar<ScalarType> data;
        std::ifstream fin;
        uint64_t stepNum;
        bool init;
    public:
        Xdatcar(std::ifstream fin_);
        Xdatcar(const Xdatcar&) = delete;
        Xdatcar(Xdatcar&&) noexcept = default;
        ~Xdatcar() = default;
        /* Operators */
        Xdatcar& operator=(Xdatcar xdatcar) noexcept;
        /* Operations */
        bool step();
        void swap(Xdatcar& __restrict xdatcar) noexcept;
        /* Getters */
        [[nodiscard]] const Poscar<ScalarType>& getCurrent() const noexcept { return data; }
        [[nodiscard]] uint64_t getStep() const noexcept { return stepNum; }
    };
}
