/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2019-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <numbers>
#include "Physica/Core/Scalar/Real.h"

namespace Physica {
    template<class> class MathConst;

    template<FloatPrec Prec>
    class MathConst<Real<Prec>> {
        using T = Real<Prec>;
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
    class PHYSICA_API MathConst<Real<FloatMP>> {
        using T = Real<FloatMP>;
        using This = MathConst<T>;
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
