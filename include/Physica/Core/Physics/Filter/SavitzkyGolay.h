/*
 * Copyright 2022-2023 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009:480-484
     */
    template<class ScalarType>
    class SavitzkyGolay {
        using MatrixType = DenseMatrix<ScalarType>;
        unsigned char lRange;
        unsigned char rRange;
        MatrixType coeffs;
        MatrixType inv;
        ScalarType delta;
    public:
        SavitzkyGolay(unsigned char lRange_, unsigned char rRange_, size_t order, ScalarType delta_);
        SavitzkyGolay(const SavitzkyGolay&) = default;
        SavitzkyGolay(SavitzkyGolay&&) noexcept = default;
        ~SavitzkyGolay() = default;
        /* Operators */
        SavitzkyGolay& operator=(SavitzkyGolay filter) noexcept { swap(filter); return *this; }
        /* Operations */
        template<class VectorType>
        void smooth(LValueVector<VectorType>& data) const;
        template<class VectorType>
        void smooth_zero(LValueVector<VectorType>& data) const;
        void swap(SavitzkyGolay& __restrict filter) noexcept;
        /* Getters */
        [[nodiscard]] unsigned char getLRange() const noexcept { return lRange; }
        [[nodiscard]] unsigned char getRRange() const noexcept { return rRange; }
        [[nodiscard]] size_t getOrder() const noexcept { return coeffs.getRow() - 1; }
        [[nodiscard]] size_t getWindowSize() const noexcept { return coeffs.getRow(); }
        [[nodiscard]] ScalarType getDelta() const noexcept { return delta; }
    };

    template<class ScalarType>
    SavitzkyGolay<ScalarType>::SavitzkyGolay(unsigned char lRange_, unsigned char rRange_, size_t order, ScalarType delta_)
            : lRange(lRange_)
            , rRange(rRange_)
            , coeffs(lRange_ + rRange_ + 1, order + 1)
            , inv(order + 1, order + 1)
            , delta(delta_) {
        for (size_t major = 0; major < coeffs.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < coeffs.getMaxMinor(); ++minor) {
                if (major == 0 && minor == 0)
                    coeffs.refFromMajorMinor(0, 0) = 1;
                const size_t row = MatrixOption::rowFromMajorMinor<MatrixType>(major, minor);
                const size_t column = MatrixOption::columnFromMajorMinor<MatrixType>(major, minor);
                const int i = int(row) - int(lRange);
                coeffs.refFromMajorMinor(major, minor) = pow(delta * i, ScalarType(column));
            }
        }

        inv = coeffs.transpose() * coeffs;
        MatrixType temp = inv.inverse();
        inv.swap(temp);
    }

    template<class ScalarType>
    template<class VectorType>
    void SavitzkyGolay<ScalarType>::smooth(LValueVector<VectorType>& data) const {
        const size_t windowSize = getWindowSize();
        Vector<ScalarType> mask(windowSize);
        for (size_t i = 0; i < windowSize; ++i)
            mask[i] = inv.row(0).asVector() * coeffs.row(i).asVector();

        Vector<ScalarType> temp = data;
        for (size_t pos = lRange; pos < data.getLength() - rRange; ++pos) {
            auto window = data.getDerived().segment(pos - lRange, pos - lRange + windowSize);
            temp[pos] = mask * window;
        }
        data.getDerived() = std::move(temp);
    }
    /**
     * Append zeros if data is undefined
     */
    template<class ScalarType>
    template<class VectorType>
    void SavitzkyGolay<ScalarType>::smooth_zero(LValueVector<VectorType>& data) const {
        const size_t windowSize = getWindowSize();
        Vector<ScalarType> mask(windowSize);
        for (size_t i = 0; i < windowSize; ++i)
            mask[i] = inv.row(0).asVector() * coeffs.row(i).asVector();

        const size_t length = data.getLength();
        Vector<ScalarType> temp = data;
        /* Append zeros to pos < 0 */ {
            const size_t rRange1 = rRange + 1;
            for (size_t pos = 0; pos < lRange; ++pos) {
                auto window = data.getDerived().segment(0, pos + rRange1);
                auto mask1 = mask.segment(windowSize - window.getLength(), windowSize);
                temp[pos] = mask1 * window;
            }
        }
        for (size_t pos = lRange; pos < length - rRange; ++pos) {
            auto window = data.getDerived().segment(pos - lRange, pos - lRange + windowSize);
            temp[pos] = mask * window;
        }
        /* Append zeros to pos > length - 1 */ {
            const size_t lRange1 = lRange + 1;
            for (size_t pos = length - lRange1; pos < length; ++pos) {
                auto window = data.getDerived().segment(pos, length);
                auto mask1 = mask.segment(0, window.getLength());
                temp[pos] = mask1 * window;
            }
        }
        data.getDerived() = std::move(temp);
    }

    template<class ScalarType>
    void SavitzkyGolay<ScalarType>::swap(SavitzkyGolay<ScalarType>& __restrict filter) noexcept {
        assert(this != &filter && "[Error]: Self swap is likely a bug");
        std::swap(lRange, filter.lRange);
        std::swap(rRange, filter.rRange);
        coeffs.swap(filter.coeffs);
        inv.swap(filter.inv);
        delta.swap(filter.delta);
    }
}
