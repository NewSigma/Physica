/*
 * Copyright 2022 WeiBo He.
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

#include "RValueVector.h"
#include "Physica/Utils/CUDA/device_obj.cuh"
#include "Physica/Utils/CUDA/DeviceProp.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<RValueVector<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        constexpr static size_t SizeAtCompile = Internal::Traits<Derived>::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Internal::Traits<Derived>::MaxSizeAtCompile;
        constexpr static bool isComplex = ScalarType::isComplex;
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
    };
}

#ifdef __CUDA_ARCH__
    #include "VectorImpl/RValueVectorImpl.cuh"
#endif
