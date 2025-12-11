/*
 * Copyright 2025 Weibo He.
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
        static_assert(T::isComplex, "[Error]: Action is complex");
    private:
        const HubbardParams<Tr>& params;

        MatrixND<T> auxField;
        DiagMatrix<Tr> matsubara;
        Trv beta;
    public:
        ActionMatrix(const HubbardParams<Tr>& params, int numFreq);
        ActionMatrix(const This&) = default;
        ActionMatrix(This&&) noexcept = default;
        ~ActionMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        template<RNG R>
        void randAuxField();
        template<RNG R>
        void randAuxField(int site);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return params; }
        [[nodiscard]] int getNumFreq() const noexcept { return auxField.getRow() / 2; }
        [[nodiscard]] int getNumSite() const noexcept { return auxField.getCol(); }
        [[nodiscard]] size_t getOrder() const noexcept { return matsubara.getRow() * getNumSite(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] auto& getAuxField() noexcept { return auxField; }
    };

    template<Scalar T>
    ActionMatrix<T>::ActionMatrix(const HubbardParams<Tr>& params, int numFreq)
            : params(params)
            , auxField(numFreq * 2, params.getNumSite())
            , matsubara(numFreq * 2) {
        assert(numFreq > 0);
        auto& diag = matsubara.diag();
        for (int k = 0; k < diag.getLength(); ++k) {
            int m = k - numFreq;
            diag[k] = Trv(2 * m + 1);
        }
        diag *= MathConst<Trv>::pi / params.getBeta();
    }

    template<Scalar T>
    void ActionMatrix<T>::assign(Matrix auto&& target) const {
        target.zeros();

        const int numSite = getNumSite();
        const Tr shiftChemMu = params.getChemMu() - params.getRepelU() * 0.5;
        kronecker(UnitMatrix<Trv>(matsubara.getRow()), params.getHoppingMatrix() - UnitMatrix<Trv>(getNumSite()) * shiftChemMu).assign_add(target);
        kronecker(matsubara, UnitMatrix<Trv>(numSite) * T(0, 1)).assign_add(target);
        for (int rowFreq = 0; rowFreq < getNumFreq() * 2; ++rowFreq) {
            int offsetR = rowFreq * numSite;
            for (int colFreq = 0; colFreq <= rowFreq; ++colFreq) {
                if (rowFreq == colFreq) {
                    target.block(offsetR, numSite, offsetR, numSite).diag() += auxField.row(0);
                    continue;
                }
                int offsetC = colFreq * numSite;
                int freq = rowFreq - colFreq;
                target.block(offsetR, numSite, offsetC, numSite).diag() = auxField.row(freq);
                target.block(offsetC, numSite, offsetR, numSite).diag() = auxField.row(freq).conjugate();
            }
        }
    }

    template<Scalar T>
    T ActionMatrix<T>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getCol());
        const int numSite = getNumSite();
        const size_t rowFreq = row / numSite;
        const size_t rowSite = row % numSite;
        const size_t colFreq = col / numSite;
        const size_t colSite = col % numSite;
        bool diagFreq = rowFreq == colFreq;
        bool diagSite = rowSite == colSite;
        if (diagFreq) {
            Tr re = 0, im = 0;
            re = params.getHoppingMatrix()(rowSite, colSite);
            if (diagSite) {
                const Tr shiftChemMu = params.getChemMu() - params.getRepelU() * 0.5;
                re += auxField(0, rowSite).real() - shiftChemMu;
                im = matsubara.diag()[rowFreq];
            }
            return T(re, im);
        }

        if (diagSite) {
            T aux;
            if (rowFreq > colFreq)
                aux = auxField(rowFreq - colFreq, rowSite);
            else
                aux = auxField(colFreq - rowFreq, rowSite).conjugate();
            return aux;
        }
        return 0;
    }

    template<Scalar T>
    template<RNG R>
    void ActionMatrix<T>::randAuxField() {
        Tr factor = sqrt(params.getRepelU() / params.getBeta());
        auxField.template random_normal<R>();
        for (int site = 0; site < getNumSite(); ++site)
            auxField(0, site).imag() = 0;
        auxField.row(0) *= factor;
        auxField.bottomRows(1) *= factor / sqrt(Trv(2));
    }

    template<Scalar T>
    template<RNG R>
    void ActionMatrix<T>::randAuxField(int site) {
        Tr factor = sqrt(params.getRepelU() / params.getBeta());
        auto aux = auxField.col(site);
        aux.template random_normal<R>();
        aux[0].imag() = 0;
        aux[0] *= factor;
        aux.tail(1) *= factor / sqrt(Trv(2));
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<ActionMatrix<T>> {
        static_assert(T::isComplex, "[Error]: Action is complex");
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}
