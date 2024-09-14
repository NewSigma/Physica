/*
 * Copyright 2019-2023 Weibo He.
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

#include <iosfwd>
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    class PHYSICA_API BasicConst {
        using ScalarType = Scalar<FloatMP>;
    public:
        double ln_2;
        double ln_10;
        double ln_2_10;
        ScalarType plotPoints;
        ScalarType expectedRelativeError;
        ScalarType stepSize;
        ScalarType R_MAX;
        ScalarType _0;
        ScalarType _1;
        ScalarType Minus_1;
        ScalarType _2;
        ScalarType Minus_2;
        ScalarType _3;
        ScalarType Minus_3;
        ScalarType _4;
        ScalarType Minus_4;
        ScalarType _10;
    public:
        BasicConst(const BasicConst&) = delete;
        BasicConst(BasicConst&&) noexcept = delete;
        ~BasicConst() = default;
        /* Operators */
        BasicConst& operator=(const BasicConst&) = delete;
        BasicConst& operator=(BasicConst&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] inline static const BasicConst& getInstance();
    private:
        BasicConst();
    };

    class PHYSICA_API MathConst {
        using ScalarType = Scalar<FloatMP>;
    public:
        ScalarType PI;
        ScalarType E;
        //Here PI_2 stands by PI / 2.
        ScalarType PI_2;
        ScalarType Minus_PI_2;
    public:
        MathConst(const MathConst&) = delete;
        MathConst(MathConst&&) noexcept = delete;
        ~MathConst() = default;
        /* Operators */
        MathConst& operator=(const MathConst&) = delete;
        MathConst& operator=(MathConst&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] inline static const MathConst& getInstance();
    private:
        MathConst();

        static ScalarType calcPI(int precision);
    };

    inline const BasicConst& BasicConst::getInstance() {
        static BasicConst basicConst{};
        return basicConst;
    }

    inline const MathConst& MathConst::getInstance() {
        static MathConst mathConst{};
        return mathConst;
    }
}
