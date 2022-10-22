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

#include "Physica/Core/Physics/ElectronicStructure/DFT/KSpaceGrid.h"
#include "FFT.h"

namespace Physica::Core {
    template<class ScalarType>
    class FFTGrid : private KSpaceGrid<ScalarType> {
        static_assert(ScalarType::isComplex == true);
    public:
        using Base = KSpaceGrid<ScalarType>;
        using RealType = typename ScalarType::RealType;
    public:
        FFTGrid() = default;
        FFTGrid(const FFT<RealType, 3>& fft);
        FFTGrid(const FFTGrid&) = default;
        FFTGrid(FFTGrid&&) noexcept = default;
        ~FFTGrid() = default;
        /* Operators */
        FFTGrid& operator=(FFTGrid grid) noexcept;
        [[nodiscard]] ScalarType& operator()(ssize_t x, ssize_t y, ssize_t z);
        [[nodiscard]] const ScalarType& operator()(ssize_t x, ssize_t y, ssize_t z) const;
        /* Operations */
        void swap(FFTGrid& grid) noexcept { Base::swap(grid); }
        /* Getters */
        using Base::asVector;
        using Base::getSize;
        using Base::getDimX;
        using Base::getDimY;
        using Base::getDimZ;
        using Base::getDim;
    };

    template<class ScalarType>
    FFTGrid<ScalarType>::FFTGrid(const FFT<RealType, 3>& fft)
            : Base(fft.getKSpace(), fft.getKSpaceSize(0) / 2, fft.getKSpaceSize(1) / 2, fft.getKSpaceSize(2) - 1) {}

    template<class ScalarType>
    FFTGrid<ScalarType>& FFTGrid<ScalarType>::operator=(FFTGrid<ScalarType> grid) noexcept {
        swap(grid);
        return *this;
    }

    template<class ScalarType>
    ScalarType& FFTGrid<ScalarType>::operator()(ssize_t x, ssize_t y, ssize_t z) {
        const ssize_t x1 = x >= 0 ? (x - getDimX()) : (x + getDimX() + 1);
        const ssize_t y1 = y >= 0 ? (y - getDimY()) : (y + getDimY() + 1);
        return Base::operator()(x1, y1, z);
    }
    
    template<class ScalarType>
    const ScalarType& FFTGrid<ScalarType>::operator()(ssize_t x, ssize_t y, ssize_t z) const {
        const ssize_t x1 = x >= 0 ? (x - getDimX()) : (x + getDimX() + 1);
        const ssize_t y1 = y >= 0 ? (y - getDimY()) : (y + getDimY() + 1);
        return Base::operator()(x1, y1, z);
    }
}
