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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h"

namespace Physica {
    /**
     * Reference:
     * [1] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016
     */
    template<Scalar T>
    class ImagKinetic {
        using This = ImagKinetic<T>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
    public:
        using GreenArray = Array2D<DenseMatrix<T>, MatrixOption::Col, 2>;
    private:
        DenseMatrix<Trv> aux;
        GreenArray greens;
        Tv sgnDet = 1;
    public:
        ImagKinetic(int numSite, int numSplit);
        ImagKinetic(const This&) = default;
        ImagKinetic(This&&) noexcept = default;
        ~ImagKinetic() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        Vector2D<Tr> calcDelta(int site, int split, Tr alpha) const noexcept;
        Vector2D<Tv> calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept;
        [[nodiscard]] Trv calcP(int site, int split, Tr alpha) const noexcept;
        void single_flip(int site, int split, Tr alpha) noexcept;

        template<RNG R>
        void random_uniform();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return aux.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return aux.getCol(); }
        [[nodiscard]] const auto& getAuxField() const noexcept { return aux; }
        [[nodiscard]] auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Tv getSign() const noexcept { return sgnDet; }
    private:
        void flipGreens(int site, int split, Vector2D<Tv> deltaRatios);
    };

    template<Scalar T>
    ImagKinetic<T>::ImagKinetic(int numSite, int numSplit)
            : aux(numSite, numSplit)
            , greens(2, numSplit, numSite, numSite) {}

    template<Scalar T>
    auto ImagKinetic<T>::calcP(int site, int split, Tr alpha) const noexcept -> Trv {
        const auto deltas = calcDelta(site, split, alpha);
        const auto ratios = calcRatio(site, split, deltas);
        return abs(ratios.prod());
    }

    template<Scalar T>
    void ImagKinetic<T>::single_flip(int site, int split, Tr alpha) noexcept {
        const auto deltas = calcDelta(site, split, alpha);
        const auto ratios = calcRatio(site, split, deltas);
        auto spins = aux.row(site);
        spins[split] = -spins[split];

        flipGreens(site, split, divide(deltas, ratios));
        Tv prod = ratios.prod();
        sgnDet *= unit(prod);
    }

    template<Scalar T>
    auto ImagKinetic<T>::calcDelta(int site, int split, Tr alpha) const noexcept -> Vector2D<Tr> {
        const Tr x = Trv(2) * alpha * aux(site, split);
        return exp(Vector2D<Tr>{-x, x}) - Trv(1); // The Eq. above Eq.(7.33) of [1]
    }

    template<Scalar T>
    auto ImagKinetic<T>::calcRatio(int site, int split, Vector2D<Tr> deltas) const noexcept -> Vector2D<Tv> {
        assert(site < getNumSite() && split < getNumSplit());
        const int split1 = (split + 1) % getNumSplit();
        Vector2D<Tv> result = deltas;
        result[0] *= Trv(1) - greens(0, split1)(site, site);
        result[1] *= Trv(1) - greens(1, split1)(site, site);
        result += Trv(1);
        return result; // Eq.(7.36) of [1]
    }

    template<Scalar T>
    void ImagKinetic<T>::flipGreens(int site, int split, Vector2D<Tv> deltaRatios) {
        // Eq. (7.44) of [1]
        const int numSite = getNumSite();
        VectorND<T> vc(numSite);
        VectorND<T> vr(numSite);
        auto flipGreen = [site, &vc, &vr](DenseMatrix<T>& green, Tv deltaRatio) {
            vc = (green - UnitMatrix<T>(green)).col(site);
            vr = green.row(site);
            green += deltaRatio * (vc * vr.transpose());
        };
        const int split1 = (split + 1) % getNumSplit();
        flipGreen(greens(0, split1), deltaRatios[0]);
        flipGreen(greens(1, split1), deltaRatios[1]);
    }

    template<Scalar T>
    template<RNG R>
    void ImagKinetic<T>::random_uniform() {
        aux.template random_uniform<R>();
        aux = unit_elem(aux - Tr(0.5));
    }

    template<Scalar T>
    void ImagKinetic<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        aux.swap(obj.aux);
        greens.swap(obj.greens);
        sgnDet.swap(obj.sgnDet);
    }
}
