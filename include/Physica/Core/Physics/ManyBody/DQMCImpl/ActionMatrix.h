/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "HubbardParams.h"

namespace Physica {
    template<Scalar T>
    class ActionMatrix : public RValueMatrix<ActionMatrix<T>> {
        using This = ActionMatrix<T>;
        using Base = RValueMatrix<This>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        static_assert(T::isComplex(), "[Error]: Action is complex");
    private:
        DiagMatrix<Tr> matsubara;
        MatrixND<T> auxField;

        const HubbardParams<Tr>& params;
    public:
        ActionMatrix(const HubbardParams<Tr>& params, int numFreq, int maxBoson) noexcept;
        ActionMatrix(const HubbardParams<Tr>& params, int numFreq) noexcept;
        ActionMatrix(const This&) = default;
        ActionMatrix(This&&) noexcept = default;
        ~ActionMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_kinetic(Matrix auto&& target) const;
        void assign_potential(Matrix auto&& target) const;
        void assign_potential(Matrix auto&& target, int site) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] This hermite(this auto&& self) noexcept;

        void flip();
        template<RNG R>
        void randAuxField();
        template<RNG R>
        [[nodiscard]] T randAuxField(int freq, int site);
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return matsubara.getOrder() * getNumSite(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] auto&& getAuxField(this auto&&) noexcept;
        [[nodiscard]] int getNumFreq() const noexcept { return matsubara.getOrder() / 2; }
        [[nodiscard]] int getNumSite() const noexcept { return auxField.getCol(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return auxField.getRow(); }
        [[nodiscard]] Tr getRepelU() const noexcept { return params.getRepelU(); }
        [[nodiscard]] Tr getBeta() const noexcept { return params.getBeta(); }
        [[nodiscard]] const auto& getParams() const noexcept { return params; }
        /* Friends */
        template<Matrix, Vector> friend class GEMV;
    private:
        ActionMatrix() = default;
    };

    template<Scalar T>
    ActionMatrix<T>::ActionMatrix(const HubbardParams<Tr>& params, int numFreq, int maxBoson) noexcept
            : matsubara(numFreq * 2)
            , auxField(maxBoson, params.getNumSite())
            , params(params) {
        assert(getNumSite() > 1 && "[Error]: Reject self hopping");
        assert(numFreq > 0);
        assert((1 <= maxBoson) && (maxBoson <= 2 * numFreq) && "[Error]: maxBoson out of range");
        auto& diag = matsubara.diag();
        for (int k = 0; k < diag.getLength(); ++k) {
            int m = k - numFreq;
            diag[k] = Trv(2 * m + 1);
        }
        diag *= MathConst<Trv>::pi;
    }

    template<Scalar T>
    ActionMatrix<T>::ActionMatrix(const HubbardParams<Tr>& params, int numFreq) noexcept : This(params, numFreq, numFreq * 2) {}

    template<Scalar T>
    void ActionMatrix<T>::assign(Matrix auto&& target) const {
        // Note that we seldom change the kinetic part, separate them to customize the potential part.
        assign_kinetic(target);
        assign_potential(target);
    }

    template<Scalar T>
    void ActionMatrix<T>::assign_kinetic(Matrix auto&& target) const {
        kronecker(IdentityMatrix<Trv>(getNumSite()) * T(0, 1), matsubara).assign(target);
        kronecker(params.getHoppingMatrix() * getBeta(), IdentityMatrix<Trv>(matsubara.getOrder())).assign_add(target);
    }

    template<Scalar T>
    void ActionMatrix<T>::assign_potential(Matrix auto&& target) const {
        const int numFreq2 = getNumFreq() * 2;
        const int numSite = getNumSite();
        for (int site = 0; site < numSite; ++site) {
            int offset = site * numFreq2;
            assign_potential(target.block(offset, numFreq2, offset, numFreq2), site);
        }
    }

    template<Scalar T>
    void ActionMatrix<T>::assign_potential(Matrix auto&& target, int site) const {
        const Tr shift = params.getBeta() * fma(params.getRepelU(), Tr(-0.5), params.getChemMu());
        const int numFreq2 = getNumFreq() * 2;
        assert(target.isSquare() && target.getRow() == numFreq2);
        for (int r = 0; r < numFreq2; ++r) {
            for (int c = 0; c < r; ++c) {
                int delta = r - c;
                if (delta < getMaxBoson()) {
                    target[r, c] = auxField[delta, site];
                    target[c, r] = auxField[delta, site].conjugate();
                }
            }
            target[r, r].real() = auxField[0, site].real() - shift;
        }
    }

    template<Scalar T>
    T ActionMatrix<T>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getCol());
        const int numFreq2 = getNumFreq() * 2;
        const size_t rowSite = row / numFreq2;
        const size_t rowFreq = row % numFreq2;
        const size_t colSite = col / numFreq2;
        const size_t colFreq = col % numFreq2;
        bool diagFreq = rowFreq == colFreq;
        bool diagSite = rowSite == colSite;
        if (diagFreq) {
            Tr re = 0, im = 0;
            re = params.getHoppingMatrix()[rowSite, colSite] * params.getBeta();
            if (diagSite) {
                const Tr shift = params.getBeta() * fma(params.getRepelU(), Tr(-0.5), params.getChemMu());
                re += auxField[0, rowSite].real() - shift;
                im = matsubara.diag()[rowFreq];
            }
            return T(re, im);
        }

        if (diagSite) {
            bool upper = rowFreq > colFreq;
            auto delta = upper ? (rowFreq - colFreq) : (colFreq - rowFreq);
            if (delta < getMaxBoson()) {
                T aux = auxField[delta, rowSite];
                return upper ? aux : aux.conjugate();
            }
        }
        return 0;
    }

    template<Scalar T>
    auto ActionMatrix<T>::hermite(this auto&& self) noexcept -> This {
        This result = std::forward<decltype(self)>(self);
        result.matsubara.diag() = -result.matsubara.diag();
        return result;
    }

    template<Scalar T>
    void ActionMatrix<T>::flip() {
        auxField = -auxField;
    }

    template<Scalar T>
    template<RNG R>
    void ActionMatrix<T>::randAuxField() {
        Tr factor = sqrt(params.getRepelU() * params.getBeta());
        auxField.template random_normal<R>();
        auxField.row(0) *= factor;
        if (getMaxBoson() > 1)
            auxField.bottomRows(1) *= factor / sqrt(Trv(2));
    }

    template<Scalar T>
    template<RNG R>
    T ActionMatrix<T>::randAuxField(int freq, int site) {
        const Tr betaU = params.getBeta() * params.getRepelU();
        const Vector3D<Tr> shifts{-betaU, 0, betaU};
        const Tr factor = sqrt(betaU);
        T result = T::template random_normal<R>() * factor;
        result.real() += shifts[std::uniform_int_distribution<int>(0, 2)(R::getInstance())];
        result.imag() += shifts[std::uniform_int_distribution<int>(0, 2)(R::getInstance())];
        if (freq != 0)
            result /= sqrt(Trv(2));
        return std::exchange(auxField[freq, site], result);
    }

    template<Scalar T>
    auto&& ActionMatrix<T>::getAuxField(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).auxField;
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<ActionMatrix<T>> {
        static_assert(T::isComplex(), "[Error]: Action is complex");
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
    };
}

#include "ActionMatrixGEMV.h"
