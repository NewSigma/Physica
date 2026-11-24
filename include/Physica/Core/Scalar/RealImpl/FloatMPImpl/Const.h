/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

namespace Physica {
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
        [[nodiscard]] static const BasicConst& getInstance() noexcept;
    private:
        BasicConst();
    };
}
