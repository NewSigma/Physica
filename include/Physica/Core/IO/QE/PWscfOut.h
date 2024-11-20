/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica::Core {
    class PHYSICA_API PWscfOut {
        using ScalarType = float64;
        constexpr static size_t DefaultBufferSize = 1024; //1024 shall be enough

        std::ifstream fin;
        VectorND<ScalarType> force;
        Array<char> buffer;
    public:
        PWscfOut(const char* path, size_t numAtom);
        PWscfOut(const PWscfOut&) = delete;
        PWscfOut(PWscfOut&&) noexcept = default;
        ~PWscfOut() = default;
        /* Operators */
        PWscfOut& operator=(PWscfOut obj) noexcept;
        /* Operations */
        [[nodiscard]] VectorND<ScalarType> makeTotalForces();
        void swap(PWscfOut& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumAtom() const noexcept { return force.getLength() / 3; }
        [[nodiscard]] const VectorND<ScalarType>& getForce() const noexcept { return force; }
    private:
        void readForce();
    };
}
