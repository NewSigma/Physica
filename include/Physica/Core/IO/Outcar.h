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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    class Outcar {
        using ScalarType = Scalar<Double, false>;

        unsigned int numAtom;
        Vector<ScalarType> force;
        ScalarType internalEnergy;
    public:
        Outcar(const char* path, unsigned int numAtom_);
        Outcar(const Outcar&) = default;
        Outcar(Outcar&&) noexcept = default;
        ~Outcar() = default;
        /* Operators */
        Outcar& operator=(Outcar& outcar) noexcept;
        /* Getters */
        [[nodiscard]] const Vector<ScalarType>& getForce() const noexcept { return force; }
        [[nodiscard]] ScalarType getInternalEnergy() const noexcept { return internalEnergy; }
        /* Helpers */
        void swap(Outcar& outcar) noexcept;
    private:
        void readForce(std::ifstream& fin, Utils::Array<char>& buffer);
        void readEnergy(std::ifstream& fin, Utils::Array<char>& buffer);
    };
}
