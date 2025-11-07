/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Utils/Container/ArrayND.h"

namespace Physica {
    template<Scalar T>
    class ProbDistribution2D {
        using This = ProbDistribution2D;
        using BucketType = ArrayND<size_t, 3>;
        using VectorType = VectorND<T>;
        using MatrixType = DenseMatrix<T>;
        using MeshType = std::pair<MatrixType, MatrixType>;

        BucketType bucket;
        VectorType seperatesX;
        VectorType seperatesY;
        T repDeltaX;
        T repDeltaY;
    public:
        ProbDistribution2D(T fromX, T toX, T fromY, T toY, size_t numBinX, size_t numBinY);
        ProbDistribution2D(const This&) = default;
        ProbDistribution2D(This&&) noexcept = default;
        ~ProbDistribution2D() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] This operator+(const This& pdf);
        void operator+=(const This& pdf);
        /* Operations */
        void sample(T x, T y);
        void clear();
        [[nodiscard]] MeshType makePosition() const;
        [[nodiscard]] MatrixType makeDistribution() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const BucketType& getBucket() const noexcept { return bucket; }
        [[nodiscard]] size_t getNumBinX() const noexcept { return bucket.dim(0); }
        [[nodiscard]] size_t getNumBinY() const noexcept { return bucket.dim(1); }
        [[nodiscard]] T getFromPointX() const noexcept { return seperatesX[0]; }
        [[nodiscard]] T getFromPointY() const noexcept { return seperatesY[0]; }
        [[nodiscard]] T getToPointX() const noexcept { return *seperatesX.crbegin(); }
        [[nodiscard]] T getToPointY() const noexcept { return *seperatesY.crbegin(); }
    private:
        size_t calcNumSample() const;
    };

    template<Scalar T>
    ProbDistribution2D<T>::ProbDistribution2D(
            T fromX, T toX, T fromY, T toY, size_t numBinX, size_t numBinY)
            : bucket({numBinX, numBinY, 1}, 0)
            , seperatesX(VectorType::linspace(fromX, toX, numBinX + 1))
            , seperatesY(VectorType::linspace(fromY, toY, numBinY + 1))
            , repDeltaX(T(numBinX) / (toX - fromX))
            , repDeltaY(T(numBinY) / (toY - fromY)) {
        assert(fromX < toX);
        assert(fromY < toY);
    }

    template<Scalar T>
    ProbDistribution2D<T> ProbDistribution2D<T>::operator+(
            const This& pdf) {
        This result(*this);
        for (size_t i = 0; i < result.getNumBinX(); ++i)
            for (size_t j = 0; j < result.getNumBinY(); ++j)
                result.bucket(i, j, 0) += pdf.bucket(i, j, 0);
        return result;
    }

    template<Scalar T>
    void ProbDistribution2D<T>::operator+=(const This& pdf) {
        for (size_t i = 0; i < getNumBinX(); ++i)
            for (size_t j = 0; j < getNumBinY(); ++j)
                bucket(i, j, 0) += pdf.bucket(i, j, 0);
    }

    template<Scalar T>
    void ProbDistribution2D<T>::sample(T x, T y) {
        const long indexX = double((x - getFromPointX()) * repDeltaX);
        const long indexY = double((y - getFromPointY()) * repDeltaY);
        if (x > getFromPointX() && 0 <= indexX && size_t(indexX) < getNumBinX())
            if (y > getFromPointY() && 0 <= indexY && size_t(indexY) < getNumBinY())
                bucket(indexX, indexY, 0) += 1;
    }

    template<Scalar T>
    void ProbDistribution2D<T>::clear() {
        for (auto& elem : bucket.asArray())
            elem = 0;
    }

    template<Scalar T>
    auto ProbDistribution2D<T>::makePosition() const -> MeshType {
        const T deltaX = (getToPointX() - getFromPointX()) / T(getNumBinX());
        const T deltaY = (getToPointY() - getFromPointY()) / T(getNumBinY());
        VectorType translateX = seperatesX.head(getNumBinX()) + (deltaX * 0.5);
        VectorType translateY = seperatesY.head(getNumBinY()) + (deltaY * 0.5);
        return MatrixType::meshgrid(translateX, translateY);
    }

    template<Scalar T>
    auto ProbDistribution2D<T>::makeDistribution() const -> MatrixType {
        const T factor = repDeltaX * repDeltaY / T(calcNumSample());
        MatrixType result(getNumBinX(), getNumBinY());
        for (size_t i = 0; i < result.getCol(); ++i) {
            auto col = result.col(i);
            for (size_t j = 0; j < result.getRow(); ++j)
                col[j] = T(bucket(i, j, 0)) * factor;
        }
        return result;
    }

    template<Scalar T>
    void ProbDistribution2D<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        bucket.swap(obj.bucket);
        seperatesX.swap(obj.seperatesX);
        seperatesY.swap(obj.seperatesY);
        repDeltaX.swap(obj.repDeltaX);
        repDeltaY.swap(obj.repDeltaY);
    }

    template<Scalar T>
    size_t ProbDistribution2D<T>::calcNumSample() const {
        size_t num = 0;
        for (size_t elem : bucket.asArray())
            num += elem;
        return num;
    }
}
