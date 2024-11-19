/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

namespace Physica::Core {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:480-484
     */
    template<Scalar T>
    class SavitzkyGolay {
        using MatrixType = DenseMatrix<T>;

        unsigned char lRange;
        unsigned char rRange;
        MatrixType coeffs;
        MatrixType inv;
        T delta;
    public:
        SavitzkyGolay(unsigned char lRange_, unsigned char rRange_, size_t order, T delta_);
        SavitzkyGolay(const SavitzkyGolay&) = default;
        SavitzkyGolay(SavitzkyGolay&&) noexcept = default;
        ~SavitzkyGolay() = default;
        /* Operators */
        SavitzkyGolay& operator=(SavitzkyGolay filter) noexcept { swap(filter); return *this; }
        /* Operations */
        template<LVector V>
        void smooth(V& data) const;
        template<LVector V>
        void smooth_zero(V& data) const;
        void swap(SavitzkyGolay& __restrict filter) noexcept;
        /* Getters */
        [[nodiscard]] unsigned char getLRange() const noexcept { return lRange; }
        [[nodiscard]] unsigned char getRRange() const noexcept { return rRange; }
        [[nodiscard]] size_t getOrder() const noexcept { return coeffs.getRow() - 1; }
        [[nodiscard]] size_t getWindowSize() const noexcept { return coeffs.getRow(); }
        [[nodiscard]] T getDelta() const noexcept { return delta; }
    };

    template<Scalar T>
    SavitzkyGolay<T>::SavitzkyGolay(unsigned char lRange_, unsigned char rRange_, size_t order, T delta_)
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
                const size_t col = MatrixOption::colFromMajorMinor<MatrixType>(major, minor);
                const int i = int(row) - int(lRange);
                coeffs.refFromMajorMinor(major, minor) = pow(delta * i, T(col));
            }
        }

        inv = coeffs.transpose() * coeffs;
        MatrixType temp = inv.inverse();
        inv.swap(temp);
    }

    template<Scalar T>
    template<LVector V>
    void SavitzkyGolay<T>::smooth(V& data) const {
        const size_t windowSize = getWindowSize();
        VectorND<T> mask(windowSize);
        for (size_t i = 0; i < windowSize; ++i)
            mask[i] = inv.row(0) * coeffs.row(i);

        VectorND<T> temp = data;
        for (size_t pos = lRange; pos < data.getLength() - rRange; ++pos) {
            auto window = data.segment(pos - lRange, pos - lRange + windowSize);
            temp[pos] = mask * window;
        }
        data = std::move(temp);
    }
    /**
     * Append zeros if data is undefined
     */
    template<Scalar T>
    template<LVector V>
    void SavitzkyGolay<T>::smooth_zero(V& data) const {
        const size_t windowSize = getWindowSize();
        VectorND<T> mask(windowSize);
        for (size_t i = 0; i < windowSize; ++i)
            mask[i] = inv.row(0) * coeffs.row(i);

        const size_t length = data.getLength();
        VectorND<T> temp = data;
        /* Append zeros to pos < 0 */ {
            const size_t rRange1 = rRange + 1;
            for (size_t pos = 0; pos < lRange; ++pos) {
                auto window = data.segment(0, pos + rRange1);
                auto mask1 = mask.segment(windowSize - window.getLength(), windowSize);
                temp[pos] = mask1 * window;
            }
        }
        for (size_t pos = lRange; pos < length - rRange; ++pos) {
            auto window = data.segment(pos - lRange, pos - lRange + windowSize);
            temp[pos] = mask * window;
        }
        /* Append zeros to pos > length - 1 */ {
            const size_t lRange1 = lRange + 1;
            for (size_t pos = length - lRange1; pos < length; ++pos) {
                auto window = data.segment(pos, length);
                auto mask1 = mask.segment(0, window.getLength());
                temp[pos] = mask1 * window;
            }
        }
        data = std::move(temp);
    }

    template<Scalar T>
    void SavitzkyGolay<T>::swap(SavitzkyGolay<T>& __restrict filter) noexcept {
        assert(this != &filter && "[Error]: Self swap is likely a bug");
        std::swap(lRange, filter.lRange);
        std::swap(rRange, filter.rRange);
        coeffs.swap(filter.coeffs);
        inv.swap(filter.inv);
        delta.swap(filter.delta);
    }
}
