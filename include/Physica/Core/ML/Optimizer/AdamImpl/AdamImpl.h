/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    template<class> class AdamImpl;

    template<Scalar T>
    class AdamBase {
    public:
        struct Args {
            T lr;
            T beta1;
            T beta2;
            T epsilon;
            T decay;
        };
    public:
        virtual void step(const Args& args, void* pTarget) = 0;
        virtual void reset(const Args& args) = 0;
    };
}

#include "AdamVector.h"
#include "AdamMatrix.h"
