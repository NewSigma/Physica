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

#include "Physica/Core/Physics/ManyBody/DQMCImpl/HubbardParams.h"

namespace Physica {
    template<Scalar T>
    class GreenSampler {
        using This = GreenSampler<T>;
    protected:
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        const HubbardParams<T>& params;
    private:
        VectorND<Tv> rsigns;
        size_t cursor = 0;
    public:
        GreenSampler(const HubbardParams<T>& params, size_t numSample);
        GreenSampler(const This&) = delete;
        GreenSampler(This&&) noexcept = delete;
        ~GreenSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] Tv calcSign() const noexcept;

        void reset() { cursor = 0; }
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return params.getNumSite(); }
        [[nodiscard]] auto getRepelU() const noexcept { return params.getRepelU(); }
        [[nodiscard]] auto getBeta() const noexcept { return params.getBeta(); }
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return params.getHoppingMatrix(); }

        [[nodiscard]] const auto& getRSigns() const noexcept { return rsigns; }
        [[nodiscard]] size_t getNumSample() const noexcept { return rsigns.getLength(); }
        [[nodiscard]] size_t getCursor() const noexcept { return cursor; }
    protected:
        void sample(T rsign) noexcept;
        [[nodiscard]] T calcDensityCorr(const MatrixND<T>& green, int siteA, int siteB) const noexcept;
        [[nodiscard]] T calcDensityCorr(const MatrixND<T>& greenA, int siteA, const MatrixND<T>& greenB, int siteB) const noexcept;
    };

    template<Scalar T>
    GreenSampler<T>::GreenSampler(const HubbardParams<T>& params, size_t numSample) : params(params), rsigns(numSample) {
        assert(numSample > 0);
    }

    template<Scalar T>
    auto GreenSampler<T>::calcSign() const noexcept -> Tv {
        return rsigns.mean();
    }

    template<Scalar T>
    void GreenSampler<T>::sample(T rsign) noexcept {
        rsigns[cursor] = rsign;
        cursor = (cursor + 1) % getNumSample();
    }
    /**
     * Reference:
     * [1] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016:198
     */
    template<Scalar T>
    T GreenSampler<T>::calcDensityCorr(const MatrixND<T>& green, int siteA, int siteB) const noexcept {
        assert(siteA < getNumSite());
        assert(siteB < getNumSite());
        T rhoA = Trv(1) - green(siteA, siteA);
        if (siteA == siteB)
            return rhoA;

        T rhoB = Trv(1) - green(siteB, siteB);
        return rhoA * rhoB - green(siteA, siteB) * green(siteB, siteA);
    }

    template<Scalar T>
    T GreenSampler<T>::calcDensityCorr(const MatrixND<T>& greenA, int siteA, const MatrixND<T>& greenB, int siteB) const noexcept {
        assert(siteA < getNumSite());
        assert(siteB < getNumSite());
        return (Trv(1) - greenA(siteA, siteA)) * (Trv(1) - greenB(siteB, siteB));
    }
}
