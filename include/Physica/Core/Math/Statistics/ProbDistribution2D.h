/*
 * Copyright 2023-2024 Weibo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Grid/GridImpl/GridStorage.h>

namespace Physica::Core {
    template<class ScalarType>
    class ProbDistribution2D {
        using This = ProbDistribution2D;
        using BucketType = GridStorage<size_t>;
        using VectorType = VectorND<ScalarType>;
        using MatrixType = DenseMatrix<ScalarType>;
        using MeshType = std::pair<MatrixType, MatrixType>;

        BucketType bucket;
        VectorType seperatesX;
        VectorType seperatesY;
        ScalarType repDeltaX;
        ScalarType repDeltaY;
    public:
        ProbDistribution2D(ScalarType fromX, ScalarType toX, ScalarType fromY, ScalarType toY, size_t numBinX, size_t numBinY);
        ProbDistribution2D(const This&) = default;
        ProbDistribution2D(This&&) noexcept = default;
        ~ProbDistribution2D() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] This operator+(const This& pdf);
        void operator+=(const This& pdf);
        /* Operations */
        void sample(ScalarType x, ScalarType y);
        void clear();
        [[nodiscard]] MeshType makePosition() const;
        [[nodiscard]] MatrixType makeDistribution() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const BucketType& getBucket() const noexcept { return bucket; }
        [[nodiscard]] size_t getNumBinX() const noexcept { return bucket.getDimX(); }
        [[nodiscard]] size_t getNumBinY() const noexcept { return bucket.getDimY(); }
        [[nodiscard]] ScalarType getFromPointX() const noexcept { return seperatesX[0]; }
        [[nodiscard]] ScalarType getFromPointY() const noexcept { return seperatesY[0]; }
        [[nodiscard]] ScalarType getToPointX() const noexcept { return *seperatesX.crbegin(); }
        [[nodiscard]] ScalarType getToPointY() const noexcept { return *seperatesY.crbegin(); }
    private:
        size_t calcNumSample() const;
    };

    template<class ScalarType>
    ProbDistribution2D<ScalarType>::ProbDistribution2D(
            ScalarType fromX, ScalarType toX, ScalarType fromY, ScalarType toY, size_t numBinX, size_t numBinY)
            : bucket({numBinX, numBinY, 1}, 0)
            , seperatesX(VectorType::linspace(fromX, toX, numBinX + 1))
            , seperatesY(VectorType::linspace(fromY, toY, numBinY + 1))
            , repDeltaX(ScalarType(numBinX) / (toX - fromX))
            , repDeltaY(ScalarType(numBinY) / (toY - fromY)) {
        assert(fromX < toX);
        assert(fromY < toY);
    }

    template<class ScalarType>
    ProbDistribution2D<ScalarType> ProbDistribution2D<ScalarType>::operator+(
            const This& pdf) {
        This result(*this);
        for (size_t i = 0; i < result.bucket.getDimX(); ++i)
            for (size_t j = 0; j < result.bucket.getDimY(); ++j)
                result.bucket(i, j, 0) += pdf.bucket(i, j, 0);
        return result;
    }

    template<class ScalarType>
    void ProbDistribution2D<ScalarType>::operator+=(const This& pdf) {
        for (size_t i = 0; i < bucket.getDimX(); ++i)
            for (size_t j = 0; j < bucket.getDimY(); ++j)
                bucket(i, j, 0) += pdf.bucket(i, j, 0);
    }

    template<class ScalarType>
    void ProbDistribution2D<ScalarType>::sample(ScalarType x, ScalarType y) {
        const long indexX = double((x - getFromPointX()) * repDeltaX);
        const long indexY = double((y - getFromPointY()) * repDeltaY);
        if (x > getFromPointX() && 0 <= indexX && size_t(indexX) < getNumBinX())
            if (y > getFromPointY() && 0 <= indexY && size_t(indexY) < getNumBinY())
                bucket(indexX, indexY, 0) += 1;
    }

    template<class ScalarType>
    void ProbDistribution2D<ScalarType>::clear() {
        for (auto& elem : bucket)
            elem = 0;
    }

    template<class ScalarType>
    typename ProbDistribution2D<ScalarType>::MeshType
    ProbDistribution2D<ScalarType>::makePosition() const {
        const ScalarType deltaX = (getToPointX() - getFromPointX()) / ScalarType(getNumBinX());
        const ScalarType deltaY = (getToPointY() - getFromPointY()) / ScalarType(getNumBinY());
        VectorType translateX = seperatesX.head(getNumBinX()) + (deltaX * 0.5);
        VectorType translateY = seperatesY.head(getNumBinY()) + (deltaY * 0.5);
        return MatrixType::meshgrid(translateX, translateY);
    }

    template<class ScalarType>
    typename ProbDistribution2D<ScalarType>::MatrixType
    ProbDistribution2D<ScalarType>::makeDistribution() const {
        const ScalarType factor = repDeltaX * repDeltaY / ScalarType(calcNumSample());
        MatrixType result(getNumBinX(), getNumBinY());
        for (size_t i = 0; i < result.getCol(); ++i) {
            auto col = result.col(i);
            for (size_t j = 0; j < result.getRow(); ++j)
                col[j] = ScalarType(bucket(i, j, 0)) * factor;
        }
        return result;
    }

    template<class ScalarType>
    void ProbDistribution2D<ScalarType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        bucket.swap(obj.bucket);
        seperatesX.swap(obj.seperatesX);
        seperatesY.swap(obj.seperatesY);
        repDeltaX.swap(obj.repDeltaX);
        repDeltaY.swap(obj.repDeltaY);
    }

    template<class ScalarType>
    size_t ProbDistribution2D<ScalarType>::calcNumSample() const {
        size_t num = 0;
        for (size_t elem : bucket)
            num += elem;
        return num;
    }
}
