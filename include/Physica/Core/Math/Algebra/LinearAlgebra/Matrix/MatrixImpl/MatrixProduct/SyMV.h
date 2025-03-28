/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

namespace Physica {
    template<Scalar T, size_t Order, Vector V>
    class SyMV : public RValueVector<SyMV<T, Order, V>> {
        using MatrixType = DenseSymmMatrix<T, Order>;
        using This = SyMV<T, Order, V>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tv;
    private:
        const MatrixType& mat;
        const V& vec;
    public:
        SyMV(const MatrixType& mat_, const V& vec_);
        SyMV(const This&) = delete;
        SyMV(This&&) noexcept = default;
        ~SyMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V1, class Executor = SeqExecutor>
        inline void assign(V1& target) const;

        [[nodiscard]] inline CoDiff<ScalarType> calc(size_t index) const;
        [[nodiscard]] inline Tv calc_value(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Scalar T, size_t Order, Vector V>
    SyMV<T, Order, V>::SyMV(const MatrixType& mat_, const V& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Scalar T, size_t Order, Vector V>
    template<Vector V1, class Executor>
    inline void SyMV<T, Order, V>::assign(V1& target) const {
        const size_t length = getLength();
        assert(length == target.getLength());
        if (length >= 16) {
            target = T(0);
            for (size_t i = 0; i < length; ++i) {
                const size_t diag = mat.toIndex1D(i, i);
                const auto seg1 = mat.asVector().segment(diag, diag + length - i);
                const auto seg2 = vec.segment(i, length);
                target[i] += seg1 * seg2;

                if (i + 1 < length) {
                    auto seg = target.segment(i + 1, length);
                    seg += seg1.tail(1) * vec.calc(i);
                }
            }
        }
        else
            Base::assign(target);
    }

    template<Scalar T, size_t Order, Vector V>
    inline auto SyMV<T, Order, V>::calc(size_t index) const -> CoDiff<ScalarType> {
        return mat.row(index) * vec;
    }

    template<Scalar T, size_t Order, Vector V>
    inline auto SyMV<T, Order, V>::calc_value(size_t index) const -> Tv {
        return mat.values().row(index) * vec.values();
    }
}

namespace Physica {
    template<Scalar T, size_t Order, Vector V>
    class Traits<SyMV<T, Order, V>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename V::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
