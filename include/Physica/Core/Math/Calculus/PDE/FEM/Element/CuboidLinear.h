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

#include "Element.h"

namespace Physica {
    template<Scalar T>
    class CuboidLinear : public Element<CuboidLinear<T>> {
        using Base = Element<CuboidLinear<T>>;
    public:
        using typename Base::VectorType;
        using typename Base::MatrixType;
        using typename Base::IndexArray;
        enum {
            LeftFrontBottom = 0,
            LeftFrontTop = 1,
            LeftBehindBottom = 2,
            LeftBehindTop = 3,
            RightFrontBottom = 4,
            RightFrontTop = 5,
            RightBehindBottom = 6,
            RightBehindTop = 7
        };
    private:
        VectorType leftFrontBottom;
        VectorType rightBehindTop;
    public:
        CuboidLinear() = default;
        CuboidLinear(VectorType leftFrontBottom_, VectorType rightBehindTop_, IndexArray globalNodes);
        CuboidLinear(const CuboidLinear&) = default;
        CuboidLinear(CuboidLinear&&) noexcept = default;
        ~CuboidLinear() = default;
        /* Operators */
        CuboidLinear& operator=(CuboidLinear obj) noexcept;
        /* Operations */
        [[nodiscard]] bool contains(const VectorType& point) const;
        [[nodiscard]] VectorType getNodePos(size_t localNode) const;
        [[nodiscard]] VectorType toLocalPos(VectorType globalPos) const;
        [[nodiscard]] VectorType toGlobalPos(VectorType localPos) const;
        void swap(CuboidLinear& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const VectorType& getLeftFrontBottom() const noexcept { return leftFrontBottom; }
        [[nodiscard]] const VectorType& getRightBehindTop() const noexcept { return rightBehindTop; }
        /* Static members */
        [[nodiscard]] static T baseFunc(size_t localNode, VectorType p);
        [[nodiscard]] static T dBase_dr(size_t localNode);
        [[nodiscard]] static T dBase_ds(size_t localNode);
        [[nodiscard]] static T dBase_dt(size_t localNode);
        [[nodiscard]] static T dBase_dr(size_t localNode, [[maybe_unused]] VectorType p) { return dBase_dr(localNode); }
        [[nodiscard]] static T dBase_ds(size_t localNode, [[maybe_unused]] VectorType p) { return dBase_ds(localNode); }
        [[nodiscard]] static T dBase_dt(size_t localNode, [[maybe_unused]] VectorType p) { return dBase_dt(localNode); }
    };

    template<Scalar T>
    CuboidLinear<T>::CuboidLinear(VectorType leftFrontBottom_, VectorType rightBehindTop_, IndexArray globalNodes)
            : Base(std::move(globalNodes))
            , leftFrontBottom(std::move(leftFrontBottom_))
            , rightBehindTop(std::move(rightBehindTop_)) {}

    template<Scalar T>
    CuboidLinear<T>& CuboidLinear<T>::operator=(CuboidLinear obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    bool CuboidLinear<T>::contains(const VectorType& point) const {
        return leftFrontBottom[0] <= point[0]
            && point[0] <= rightBehindTop[0]
            && leftFrontBottom[1] <= point[1]
            && point[1] <= rightBehindTop[1]
            && leftFrontBottom[2] <= point[2]
            && point[2] <= rightBehindTop[2];
    }

    template<Scalar T>
    CuboidLinear<T>::VectorType CuboidLinear<T>::getNodePos(size_t localNode) const {
        switch (localNode) {
            case LeftFrontBottom:
                return leftFrontBottom;
            case LeftFrontTop:
                return {leftFrontBottom[0], leftFrontBottom[1], rightBehindTop[2]};
            case LeftBehindBottom:
                return {leftFrontBottom[0], rightBehindTop[1], leftFrontBottom[2]};
            case LeftBehindTop:
                return {leftFrontBottom[0], rightBehindTop[1], rightBehindTop[2]};
            case RightFrontBottom:
                return {rightBehindTop[0], leftFrontBottom[1], leftFrontBottom[2]};
            case RightFrontTop:
                return {rightBehindTop[0], leftFrontBottom[1], rightBehindTop[2]};
            case RightBehindBottom:
                return {rightBehindTop[0], rightBehindTop[1], leftFrontBottom[2]};
            case RightBehindTop:
                return rightBehindTop;
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<Scalar T>
    CuboidLinear<T>::VectorType CuboidLinear<T>::toLocalPos(VectorType globalPos) const {
        return divide(T(2) * globalPos - rightBehindTop - leftFrontBottom, rightBehindTop - leftFrontBottom);
    }

    template<Scalar T>
    CuboidLinear<T>::VectorType CuboidLinear<T>::toGlobalPos(VectorType localPos) const {
        return (leftFrontBottom + rightBehindTop - hadamard(rightBehindTop - leftFrontBottom, localPos)) * T(0.5);
    }

    template<Scalar T>
    void CuboidLinear<T>::swap(CuboidLinear& __restrict obj) noexcept {
        Base::swap(obj);
        leftFrontBottom.swap(obj.leftFrontBottom);
        rightBehindTop.swap(obj.rightBehindTop);
    }

    template<Scalar T>
    T CuboidLinear<T>::baseFunc(size_t localNode, VectorType p) {
        switch (localNode) {
            case LeftFrontBottom:
                return (T(1) - p[0]) * (T(1) - p[1]) * (T(1) - p[2]) * T(0.125);
            case LeftFrontTop:
                return (T(1) - p[0]) * (T(1) - p[1]) * (T(1) + p[2]) * T(0.125);
            case LeftBehindBottom:
                return (T(1) - p[0]) * (T(1) + p[1]) * (T(1) - p[2]) * T(0.125);
            case LeftBehindTop:
                return (T(1) - p[0]) * (T(1) + p[1]) * (T(1) + p[2]) * T(0.125);
            case RightFrontBottom:
                return (T(1) + p[0]) * (T(1) - p[1]) * (T(1) - p[2]) * T(0.125);
            case RightFrontTop:
                return (T(1) + p[0]) * (T(1) - p[1]) * (T(1) + p[2]) * T(0.125);
            case RightBehindBottom:
                return (T(1) + p[0]) * (T(1) + p[1]) * (T(1) - p[2]) * T(0.125);
            case RightBehindTop:
                return (T(1) + p[0]) * (T(1) + p[1]) * (T(1) + p[2]) * T(0.125);
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<Scalar T>
    T CuboidLinear<T>::dBase_dr(size_t localNode) {
        switch (localNode) {
            case LeftFrontBottom:
            case LeftFrontTop:
            case LeftBehindBottom:
            case LeftBehindTop:
                return T(-0.125);
            case RightFrontBottom:
            case RightFrontTop:
            case RightBehindBottom:
            case RightBehindTop:
                return T(0.125);
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<Scalar T>
    T CuboidLinear<T>::dBase_ds(size_t localNode) {
        switch (localNode) {
            case LeftFrontBottom:
            case LeftFrontTop:
            case RightFrontBottom:
            case RightFrontTop:
                return T(-0.125);
            case LeftBehindBottom:
            case LeftBehindTop:
            case RightBehindBottom:
            case RightBehindTop:
                return T(0.125);
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<Scalar T>
    T CuboidLinear<T>::dBase_dt(size_t localNode) {
        switch (localNode) {
            case LeftFrontBottom:
            case LeftBehindBottom:
            case RightFrontBottom:
            case RightBehindBottom:
                return T(-0.125);
            case LeftFrontTop:
            case LeftBehindTop:
            case RightFrontTop:
            case RightBehindTop:
                return T(0.125);
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<CuboidLinear<T>> {
    public:
        constexpr static unsigned int Dim = 3;
        constexpr static unsigned int Order = 1;
        constexpr static unsigned int NumPoint = 8;
        constexpr static unsigned int DegreeOfFreedom = NumPoint * Order;
        using ScalarType = T;
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Dim, Dim>;
    };
}
