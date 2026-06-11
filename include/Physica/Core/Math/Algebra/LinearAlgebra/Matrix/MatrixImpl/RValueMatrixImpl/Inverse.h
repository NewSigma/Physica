/*
 * Copyright 2021-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar, bool> class DenseLU;

    template<Matrix M>
    class Inverse<M> : public RValueMatrix<Inverse<M>> {
        using This = Inverse<M>;
        using Base = RValueMatrix<This>;
        using typename Base::T;

        decay_rvalue_t<M> mat;
    public:
        explicit Inverse(M&& mat_);
        Inverse(const This&) = default;
        Inverse(This&&) noexcept = default;
        ~Inverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const M& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    private:
        void assign2D(Matrix auto& target) const;
        void assign3D(Matrix auto& target) const;
    };

    template<Matrix M>
    Inverse<M>::Inverse(M&& mat_) : mat(std::forward<M>(mat_)) {
        assert(mat.getRow() == mat.getCol());
    }

    template<Matrix M>
    void Inverse<M>::assign(Matrix auto& target) const {
        using M2 = std::remove_cvref_t<decltype(target)>;
        constexpr size_t Order = std::max(Base::getRowAtCompile(), M2::getRowAtCompile());
        if constexpr (Order == 1)
            target = reciprocal(mat[0, 0]);
        else if constexpr (Order == 2)
            assign2D(target);
        else if constexpr (Order == 3)
            assign3D(target);
        else {
            auto lu = DenseLU<T, true>(mat);
            target = lu.getMatrixL().inv();
            target *= lu.getPerm().inv();
            target = lu.getMatrixU().inv() * target;
        }
    }

    template<Matrix M>
    void Inverse<M>::assign2D(Matrix auto& target) const {
        target[0, 0] = mat[1, 1];
        target[0, 1] = -mat[0, 1];
        target[1, 0] = -mat[1, 0];
        target[1, 1] = mat[0, 0];
        target *= reciprocal(mat.det());
    }

    template<Matrix M>
    void Inverse<M>::assign3D(Matrix auto& target) const {
        if constexpr (MatrixMajor::isRowMatrix<decltype(target)>()) {
            target[0, 0] = (mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1]);
            target[0, 1] = (mat[2, 1] * mat[0, 2] - mat[0, 1] * mat[2, 2]);
            target[0, 2] = (mat[0, 1] * mat[1, 2] - mat[1, 1] * mat[0, 2]);
            target[1, 0] = (mat[2, 0] * mat[1, 2] - mat[1, 0] * mat[2, 2]);
            target[1, 1] = (mat[0, 0] * mat[2, 2] - mat[2, 0] * mat[0, 2]);
            target[1, 2] = (mat[1, 0] * mat[0, 2] - mat[0, 0] * mat[1, 2]);
            target[2, 0] = (mat[1, 0] * mat[2, 1] - mat[2, 0] * mat[1, 1]);
            target[2, 1] = (mat[2, 0] * mat[0, 1] - mat[0, 0] * mat[2, 1]);
            target[2, 2] = (mat[0, 0] * mat[1, 1] - mat[1, 0] * mat[0, 1]);
        }
        else {
            target[0, 0] = (mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1]);
            target[1, 0] = (mat[2, 0] * mat[1, 2] - mat[1, 0] * mat[2, 2]);
            target[2, 0] = (mat[1, 0] * mat[2, 1] - mat[2, 0] * mat[1, 1]);
            target[0, 1] = (mat[2, 1] * mat[0, 2] - mat[0, 1] * mat[2, 2]);
            target[1, 1] = (mat[0, 0] * mat[2, 2] - mat[2, 0] * mat[0, 2]);
            target[2, 1] = (mat[2, 0] * mat[0, 1] - mat[0, 0] * mat[2, 1]);
            target[0, 2] = (mat[0, 1] * mat[1, 2] - mat[1, 1] * mat[0, 2]);
            target[1, 2] = (mat[1, 0] * mat[0, 2] - mat[0, 0] * mat[1, 2]);
            target[2, 2] = (mat[0, 0] * mat[1, 1] - mat[1, 0] * mat[0, 1]);
        }
        target *= reciprocal(mat.det());
    }

    template<Matrix M>
    __host__ __device__ consteval int Inverse<M>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Inverse<M>> : public Traits<M> {
    public:
        using ExprType = M;
    };
}
