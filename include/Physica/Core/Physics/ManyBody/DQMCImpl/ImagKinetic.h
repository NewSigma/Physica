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
        using GreenPair = Array<MatrixND<T>, 2>;
    private:
        MatrixND<Trv> aux;
        GreenPair greens;
        Tv rsign = 1;
    public:
        ImagKinetic(int numSite, int numSplit);
        ImagKinetic(const This&) = default;
        ImagKinetic(This&&) noexcept = default;
        ~ImagKinetic() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        Vector2D<Tr> calcDelta(int site, int split, Tr alpha) const noexcept;
        Vector2D<Tv> calcRatio(int site, Vector2D<Tr> deltas) const noexcept;
        [[nodiscard]] Trv calcP(int site, int split, Tr alpha) const noexcept;
        void single_flip(int site, int split) noexcept;
        void single_flip(int site, int split, Tr alpha) noexcept;

        template<RNG R>
        void random_uniform();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getAuxField() const noexcept { return aux; }
        [[nodiscard]] int getNumSite() const noexcept { return aux.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return aux.getCol(); }
        [[nodiscard]] auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Tv getRSign() const noexcept;
    private:
        void flipGreens(int site, Vector2D<Tv> deltaRatios);
    };

    template<Scalar T>
    ImagKinetic<T>::ImagKinetic(int numSite, int numSplit)
            : aux(numSite, numSplit)
            , greens(2, numSite, numSite) {}

    template<Scalar T>
    auto ImagKinetic<T>::calcDelta(int site, int split, Tr alpha) const noexcept -> Vector2D<Tr> {
        const Tr x = Trv(2) * alpha * aux(site, split);
        return exp(Vector2D<Tr>{-x, x}) - Trv(1); // The Eq. above Eq.(7.33) of [1]
    }

    template<Scalar T>
    auto ImagKinetic<T>::calcRatio(int site, Vector2D<Tr> deltas) const noexcept -> Vector2D<Tv> {
        assert(site < getNumSite());
        Vector2D<Tv> result = deltas;
        for (int spin : {0, 1})
            result[spin] *= Trv(1) - greens[spin](site, site);
        result += Trv(1);
        return result; // Eq.(7.36) of [1]
    }

    template<Scalar T>
    auto ImagKinetic<T>::calcP(int site, int split, Tr alpha) const noexcept -> Trv {
        const auto deltas = calcDelta(site, split, alpha);
        const auto ratios = calcRatio(site, deltas);
        return abs(ratios.prod());
    }

    template<Scalar T>
    void ImagKinetic<T>::single_flip(int site, int split) noexcept {
        auto spins = aux.row(site);
        spins[split] = -spins[split];
    }

    template<Scalar T>
    void ImagKinetic<T>::single_flip(int site, int split, Tr alpha) noexcept {
        const auto deltas = calcDelta(site, split, alpha);
        const auto ratios = calcRatio(site, deltas);
        single_flip(site, split);

        flipGreens(site, divide(deltas, ratios));
        rsign *= unit(ratios.prod());
    }
    /**
     * Note that observable is <A> = <As>/<s> and signs may differ a factor. We are interested in relative sign only.
     */
    template<Scalar T>
    auto ImagKinetic<T>::getRSign() const noexcept -> Tv {
        return rsign;
    }

    template<Scalar T>
    void ImagKinetic<T>::flipGreens(int site, Vector2D<Tv> deltaRatios) {
        // Eq. (7.44) of [1]
        const int numSite = getNumSite();
        VectorND<T> vc(numSite);
        VectorND<T> vr(numSite);
        for (int spin : {0, 1}) {
            auto& green = greens[spin];
            vc = (green - UnitMatrix<T>(green)).col(site);
            vr = green.row(site);
            green += deltaRatios[spin] * (vc * vr.transpose());
        }
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
        rsign.swap(obj.rsign);
    }
}
