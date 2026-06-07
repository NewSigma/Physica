/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Random.h"

namespace Physica {
    /**
    * Reference:
    * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:406-408
    */
    class PHYSICA_API Sobol {
        using This = Sobol;
        constexpr static int8_t Degree[]{1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5};
        constexpr static int16_t PolyCoeff[]{0, 1, 1, 2, 1, 4, 2, 4, 7, 11, 13, 14};
        constexpr static int MaxBit = 30, MaxDim = sizeof(Degree) / sizeof(int8_t);
        static_assert(MaxDim == sizeof(PolyCoeff) / sizeof(int16_t), "[Error]: Data is not complete");

        unsigned int numStep;
        Array<unsigned int> mask;
        Array<unsigned int> buffer;
    public:
        Sobol();
        Sobol(const This&) = default;
        Sobol(This&&) noexcept = default;
        ~Sobol() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(int i);
        template<Vector V>
        void fill(V& x);
        void reset();

        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
    public:
        int pre_step();
    };

    template<Vector V>
    void Sobol::fill(V& x) {
        static_assert(x.getSizeAtCompile() <= MaxDim, "[Error]: Vector is too long");
        constexpr static typename V::ScalarType factor = 1.0 / (1UL << MaxBit);
        assert(x.getLength() <= MaxDim && "[Error]: Vector is too long");
        const int m = pre_step();
        for (int i = 0; i < x.getLength(); ++i) {
            buffer[i] ^= mask[m + i];
            x[i] = buffer[i];
        }
        x *= factor;
    }
}
