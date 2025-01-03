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
#include "Physica/Core/Scalar/Real.h"

namespace Physica::Core {
    /*!
     * Basic consts that initialize directly.
     */
    BasicConst::BasicConst()
            : ln_2(std::log(2))
            , ln_10(std::log(10))
            , ln_2_10(std::log(2) / std::log(10))
            , plotPoints(static_cast<SignedMPUnit>(20)) {
        expectedRelativeError = T(1, 1 - GlobalPrecision);
        expectedRelativeError.setByte(0, 1);

        stepSize = T(1, - GlobalPrecision / 2); //(- GlobalPrecision / 2) is selected by experience.
        stepSize.setByte(0, 1);

        R_MAX = 2147483647;
        _0 = 0;
        _1 = 1;
        Minus_1 = -1;
        _2 = 2;
        Minus_2 = -2;
        _3 = 3;
        Minus_3 = -3;
        _4 = 4;
        Minus_4 = -4;
        _10 = 10;
    }

    const BasicConst& BasicConst::getInstance() noexcept {
        static BasicConst basicConst{};
        return basicConst;
    }
}
