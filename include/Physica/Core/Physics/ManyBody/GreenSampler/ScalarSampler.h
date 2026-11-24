/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Scalar/Diff.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "GreenSampler.h"
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ImagKinetic.h"

namespace Physica {
    template<Scalar T>
    class ScalarSampler : public GreenSampler<T> {
        using This = ScalarSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;

        using typename Base::Tv;
        using typename Base::Trv;
    public:
        enum Observable : char {
            Density,
            Density2,
            DoubleOccupy,
            MagMoment,
            Kinetic,
            Potential,
            Internal,
        };
    private:
        VectorND<T> observes;
    public:
        ScalarSampler(const HubbardParams<T>& params, size_t numSample);
        ScalarSampler(const This&) = delete;
        ScalarSampler(This&&) noexcept = delete;
        ~ScalarSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(T observe, T rsign);
        void sample(const GreenPair& greens, T rsign, Observable type);

        [[nodiscard]] T calcRawMean() const;
        [[nodiscard]] T calcMean() const;
        /* Getters */
        using Base::getNumSite;
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    private:
        [[nodiscard]] T calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable typeX) const noexcept;
    };

    template<Scalar T>
    ScalarSampler<T>::ScalarSampler(const HubbardParams<T>& params, size_t numSample)
            : Base(params, numSample)
            , observes(numSample) {}

    template<Scalar T>
    void ScalarSampler<T>::sample(T observe, T rsign) {
        observes[Base::getCursor()] = observe * rsign;
        Base::sample(rsign);
    }

    template<Scalar T>
    void ScalarSampler<T>::sample(const GreenPair& greens, T rsign, Observable type) {
        sample(calcObservable(greens[0], greens[1], type), rsign);
    }

    template<Scalar T>
    T ScalarSampler<T>::calcRawMean() const {
        if constexpr (T::isDiffable())
            return Base::calcGradScore(observes.mean());
        else
            return observes.mean();
    }

    template<Scalar T>
    T ScalarSampler<T>::calcMean() const {
        return calcRawMean() / Base::calcRSign();
    }

    template<Scalar T>
    T ScalarSampler<T>::calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable typeX) const noexcept {
        switch (typeX) {
        case Density:
            return T(2) - greenU.diag().reals().mean() - greenD.diag().reals().mean();
        case Density2: {
            int numSite = getNumSite();
            T result = 0;
            for (int siteA = 0; siteA < numSite; ++siteA) {
                for (int siteB = 0; siteB < numSite; ++siteB) {
                    result += Base::calcDensityCorr(greenU, siteA, siteB);
                    result += Base::calcDensityCorr(greenD, siteA, siteB);
                }
            }
            result /= Trv(numSite * numSite);

            T rhoU = T(1) - greenU.diag().reals().mean();
            T rhoD = T(1) - greenD.diag().reals().mean();
            result += rhoU * rhoD * Trv(2);
            return result;
        }
        case DoubleOccupy:
            return (T(1) - greenU.diag().reals()) * (T(1) - greenD.diag().reals()) / T(getNumSite());
        case MagMoment:
            return calcObservable(greenU, greenD, Density) - calcObservable(greenU, greenD, DoubleOccupy) * Trv(2);
        case Kinetic: {
            const auto& hoppingT = Base::getHoppingMatrix();
            return -(hoppingT * (greenU + greenD)).trace().real() / T(getNumSite());
        }
        case Potential:
            return calcObservable(greenU, greenD, DoubleOccupy) * Base::getRepelU();
        case Internal:
            return calcObservable(greenU, greenD, Kinetic) + calcObservable(greenU, greenD, Potential);
        default:
            unreachable();
        }
    }
}
