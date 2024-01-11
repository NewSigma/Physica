/*
 * Copyright 2019-2023 WeiBo He.
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
    //!GlobalPrecision is the max length of Scalar<MultiPrecision>
    constexpr int GlobalPrecision = 4;
    static_assert(GlobalPrecision > 1, "GlobalPrecision must be larger than 1.");

    //!RelativeError is the stop criteria of iterate method that uses Scalar<Float> or Scalar<Double>.
    constexpr double RelativeError = 1e-5;

    class PHYSICA_API BasicConst {
    public:
        double ln_2;
        double ln_10;
        double ln_2_10;
        Scalar<MultiPrecision> plotPoints;
        Scalar<MultiPrecision> expectedRelativeError;
        Scalar<MultiPrecision> stepSize;
        Scalar<MultiPrecision> R_MAX;
        Scalar<MultiPrecision> _0;
        Scalar<MultiPrecision> _1;
        Scalar<MultiPrecision> Minus_1;
        Scalar<MultiPrecision> _2;
        Scalar<MultiPrecision> Minus_2;
        Scalar<MultiPrecision> _3;
        Scalar<MultiPrecision> Minus_3;
        Scalar<MultiPrecision> _4;
        Scalar<MultiPrecision> Minus_4;
        Scalar<MultiPrecision> _10;
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
    public:
        MultiScalar PI;
        MultiScalar E;
        //Here PI_2 stands by PI / 2.
        MultiScalar PI_2;
        MultiScalar Minus_PI_2;
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

        static MultiScalar calcPI(int precision);
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
