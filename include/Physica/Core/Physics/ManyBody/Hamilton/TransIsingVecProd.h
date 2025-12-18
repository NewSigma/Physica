/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Physics/ManyBody/ReprSpace/PauliMatrix.h"
#include "TransIsingMatrix.h"

namespace Physica {
    template<Scalar T, int Dim, int NumSite, BoundaryCond BC, Vector V>
    class TransIsingVecProd : public RValueVector<TransIsingVecProd<T, Dim, NumSite, BC, V>> {
        using MatrixType = TransIsingMatrix<T, Dim, NumSite, BC>;
        using This = TransIsingVecProd<T, Dim, NumSite, BC, V>;
        using Base = RValueVector<This>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        const MatrixType& mat;
        const V& vec;
    public:
        TransIsingVecProd(const MatrixType& mat_, const V& vec_);
        TransIsingVecProd(const This&) = default;
        TransIsingVecProd(This&&) noexcept = default;
        ~TransIsingVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const V& getRHS() const noexcept { return vec; }
    private:
        /* Getters */
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
        [[nodiscard]] T getCouplingJ() const noexcept { return mat.getCouplingJ(); }
        [[nodiscard]] T getTransG() const noexcept { return mat.getTransG(); }
    };

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC, Vector V>
    TransIsingVecProd<T, Dim, NumSite, BC, V>::TransIsingVecProd(const MatrixType& mat_, const V& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC, Vector V>
    template<ExecutePolicy P>
    void TransIsingVecProd<T, Dim, NumSite, BC, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        target.zeros();

        const auto& repr = getRepr();
        for (size_t i = 0; i < getLength(); ++i) {
            using enum PauliIndex;
            auto psi0 = repr[i];
            const T g = getTransG() * vec.calc(i);
            for (unsigned int site = 0; site < NumSite; ++site) {
                auto psi = psi0;
                target[repr[psi]] -= PauliMatrix<T, X>(site).apply(psi) * g;
            }

            const T couplingJ = getCouplingJ() * vec.calc(i);
            if constexpr (BC == BoundaryCond::OBC) {
                for (unsigned int site = 0; site < NumSite - 1; ++site) {
                    auto coeff = PauliMatrix<T, Z>(site).apply(psi0);
                    coeff *= PauliMatrix<T, Z>(site + 1).apply(psi0);
                    target[i] -= couplingJ * coeff;
                }
            }
            else {
                static_assert(BC == BoundaryCond::PBC, "[Error]: Unsupported BoundaryCond");
                for (unsigned int site = 0; site < NumSite; ++site) {
                    auto coeff = PauliMatrix<T, Z>(site).apply(psi0);
                    coeff *= PauliMatrix<T, Z>((site + 1) % NumSite).apply(psi0);
                    target[i] -= couplingJ * coeff;
                }
            }
        }
    }
}

namespace Physica {
    template<Scalar T, int Dim, int NumSite, BoundaryCond BC, Vector V>
    class Traits<TransIsingVecProd<T, Dim, NumSite, BC, V>> {
        using MatrixType = TransIsingMatrix<T, Dim, NumSite, BC>;
        using T1 = V::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T, T1>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
