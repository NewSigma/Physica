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

namespace Physica::Core {
    class PHYSICA_API BasicConst {
        using T = Real<FloatMP>;
    public:
        double ln_2;
        double ln_10;
        double ln_2_10;
        T plotPoints;
        T expectedRelativeError;
        T stepSize;
        T R_MAX;
        T _0;
        T _1;
        T Minus_1;
        T _2;
        T Minus_2;
        T _3;
        T Minus_3;
        T _4;
        T Minus_4;
        T _10;
    public:
        BasicConst(const BasicConst&) = delete;
        BasicConst(BasicConst&&) noexcept = delete;
        ~BasicConst() = default;
        /* Operators */
        BasicConst& operator=(const BasicConst&) = delete;
        BasicConst& operator=(BasicConst&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] static const BasicConst& getInstance();
    private:
        BasicConst();
    };
}
