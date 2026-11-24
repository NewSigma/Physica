/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:480-484
     */
    template<Scalar T>
    class SavitzkyGolay {
        using This = SavitzkyGolay<T>;
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
        void smooth(Vector auto& data) const;
        void smooth_zero(Vector auto& data) const;
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
        const size_t maxMinor = coeffs.getMaxMinor();
        auto flat = coeffs.flatten();
        for (size_t index = 0; index < flat.getLength(); ++index) {
            const size_t major = index / maxMinor;
            const size_t minor = index % maxMinor;
            const size_t row = MatrixMajor::rowFromMajorMinor<MatrixType>(major, minor);
            const size_t col = MatrixMajor::colFromMajorMinor<MatrixType>(major, minor);
            flat[index] = pow0(delta * (int(row) - int(lRange)), T(col));
        }

        inv = coeffs.transpose() * coeffs;
        MatrixType temp = inv.inv();
        inv.swap(temp);
    }

    template<Scalar T>
    void SavitzkyGolay<T>::smooth(Vector auto& data) const {
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
    void SavitzkyGolay<T>::smooth_zero(Vector auto& data) const {
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
    void SavitzkyGolay<T>::swap(This& __restrict filter) noexcept {
        assert(this != &filter && "[Error]: Self swap is likely a bug");
        std::swap(lRange, filter.lRange);
        std::swap(rRange, filter.rRange);
        coeffs.swap(filter.coeffs);
        inv.swap(filter.inv);
        delta.swap(filter.delta);
    }
}
