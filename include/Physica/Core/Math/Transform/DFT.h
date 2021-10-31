/*
 * Copyright 2020-2021 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"

namespace Physica::Core {
    template<class ScalarType>
    class DFT {
        using RealType = typename ScalarType::ScalarType;
        using ComplexType = ComplexScalar<RealType>;
        Vector<ComplexType> data;
        ScalarType distance;
    public:
        DFT(const Vector<ScalarType>& data_, const ScalarType& distance_);
        DFT(const DFT& dft);
        DFT(DFT&& dft) noexcept;
        ~DFT() = default;
        /* Operators */
        DFT& operator=(DFT dft);
        ComplexType operator()(size_t i) { return data[i]; }
        /* Transforms */
        inline void transform();
        inline void invTransform();
        /* Getters */
        [[nodiscard]] ComplexType getComponent(ssize_t index) const;
        [[nodiscard]] Vector<ComplexType> getComponents() const;
        [[nodiscard]] const Vector<ComplexType>& getData() const noexcept { return data; }
        /* Helpers */
        void swap(DFT& dft);
    private:
        void transformImpl(const RealType& phase);
    };
}

#include "DFTImpl.h"
