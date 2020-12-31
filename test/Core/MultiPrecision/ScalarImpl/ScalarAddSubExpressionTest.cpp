/*
 * Copyright 2020 WeiBo He.
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
#include "Physica/Core/MultiPrecision/Scalar.h"

using namespace Physica::Core;

int main() {
    Scalar<MultiPrecision, false> s1(1);
    Scalar<MultiPrecision, false> s2(1);
    s2.setPower(-GlobalPrecision - 1);
    Scalar<MultiPrecision, false> s3(1);
    Scalar<MultiPrecision, false> result1(s1 + s2 - s3);
    Scalar<MultiPrecision, false> temp(s1 + s2);
    Scalar<MultiPrecision, false> result2(temp - s3);
    return result1.isZero() || !result2.isZero();
}
