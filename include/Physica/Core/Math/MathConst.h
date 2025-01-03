/*
 * Copyright 2019-2024 Weibo He.
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

#include <numbers>
#include "Physica/Core/Scalar/Real.h"

namespace Physica::Core {
    template<ScalarOption Option>
    class MathConst {
        using T = Real<Option>;
        using MachineType = T::MachineType;
    public:
        constexpr static T e = std::numbers::e_v<MachineType>;
        constexpr static T log2e = std::numbers::log2e_v<MachineType>;
        constexpr static T log10e = std::numbers::log10e_v<MachineType>;
        constexpr static T pi = std::numbers::pi_v<MachineType>;
        constexpr static T inv_pi = std::numbers::inv_pi_v<MachineType>;
        constexpr static T inv_sqrtpi = std::numbers::inv_sqrtpi_v<MachineType>;
        constexpr static T ln2 = std::numbers::ln2_v<MachineType>;
        constexpr static T ln10 = std::numbers::ln10_v<MachineType>;
        constexpr static T sqrt2 = std::numbers::sqrt2_v<MachineType>;
        constexpr static T sqrt3 = std::numbers::sqrt3_v<MachineType>;
        constexpr static T inv_sqrt3 = std::numbers::inv_sqrt3_v<MachineType>;
        constexpr static T egamma = std::numbers::egamma_v<MachineType>;
        constexpr static T phi = std::numbers::phi_v<MachineType>;
    };

    template<>
    class PHYSICA_API MathConst<FloatMP> {
        using T = Real<FloatMP>;
    public:
        T PI;
        T E;
        // Here PI_2 stands by PI / 2.
        T PI_2;
        T Minus_PI_2;
    public:
        MathConst(const MathConst&) = delete;
        MathConst(MathConst&&) noexcept = delete;
        ~MathConst() = default;
        /* Operators */
        MathConst& operator=(const MathConst&) = delete;
        MathConst& operator=(MathConst&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] static const MathConst& getInstance() noexcept;
    private:
        MathConst();

        static T calcPI(int precision);
    };
}
