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

#include <cmath>
#include <ostream>
#include <random>
#ifdef PHYSICA_CUDA
    #include <cuda/std/limits>
#endif
#include "Physica/Core/Exception/NoImplException.h"
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Scalar.h"
#include "RealImpl/ScalarBase.h"

namespace Physica::Core {
    template<Scalar, class RandomType> class GaussRandomPool;

    inline int countLeadingZeros(MPUnit n) noexcept;
    inline int countBackZeros(unsigned long n) noexcept;
    inline int countOnes(MPUnit n) noexcept;
    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte);
}

#include "Rational.h"
#include "RealImpl/RealImpl.h"
#include "RealImpl/Float32.h"
#include "RealImpl/Float64.h"
#include "RealImpl/Math.h"
#ifdef PHYSICA_CUDA
    #include "RealImpl/Float16.h"
#endif
#include "RealImpl/FloatMP.h"
#include "RealImpl/SIMD.h"
