/*
 * Copyright 2026 Weibo He.
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

#include "GPModel.h"
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"

namespace Physica {
    /**
     * \class BayesOpt handles: 1. constraints; 2. acquisition function;
     */
    template<Scalar T>
    class BayesOpt {
        using This = BayesOpt<T>;
        using Tv = T::ValueType;
        using GlobalOpt = Vegas<Tv, true>;

        constexpr static Tv Lowest = std::numeric_limits<Tv>::lowest();
    public:
        using Cube = GlobalOpt::Cube;
    private:
        GlobalOpt vegas;
        VectorND<T> optimal;
        Tv maximum;
    public:
        explicit BayesOpt(GlobalOpt vegas, VectorND<T> optimal, Tv maximum);
        BayesOpt(const This&) = default;
        BayesOpt(This&&) noexcept = default;
        ~BayesOpt() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<RNG R>
        [[nodiscard]] auto propose(std::invocable<VectorND<T>> auto fn, const instanceof<GPModel> auto& model, int iteration);

        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumFeature() const noexcept { return vegas.getDim(); }
        [[nodiscard]] const auto& getOptimal() const noexcept { return optimal; }
        [[nodiscard]] Tv getMaximum() const noexcept { return maximum; }
    private:
        [[nodiscard]] Tv calcEI(const instanceof<GPModel> auto& model, const VectorND<Tv>& x) const noexcept;
    };

    template<Scalar T>
    BayesOpt<T>::BayesOpt(GlobalOpt vegas, VectorND<T> optimal, Tv maximum)
            : vegas(std::move(vegas))
            , optimal(std::move(optimal))
            , maximum(maximum) {}

    template<Scalar T>
    template<RNG R>
    auto BayesOpt<T>::propose(std::invocable<VectorND<T>> auto fn, const instanceof<GPModel> auto& model, int iteration) {
        for (size_t d = 0; d < vegas.getDim(); ++d)
            vegas.mesh_uniform(d);

        auto xnew = vegas.template maximize<R>([this, &model](const VectorND<Tv>& x) {
            return calcEI(model, x);
        }, iteration);

        Tv ynew = fn(xnew);
        if (ynew > maximum) {
            maximum = ynew;
            xnew.assign(optimal);
        }
        return std::make_pair(std::move(xnew), std::move(ynew));
    }

    template<Scalar T>
    void BayesOpt<T>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        vegas.swap(obj.vegas);
        optimal.swap(obj.optimal);
        maximum.swap(obj.maximum);
    }

    template<Scalar T>
    auto BayesOpt<T>::calcEI(const instanceof<GPModel> auto& model, const VectorND<Tv>& x) const noexcept -> Tv {
        const auto [mean, devia] = model.predict(x);
        if (devia.isSubNormal())
            return 0;
        Tv zscore = (mean - maximum) / devia;
        // Coeff has been ignored without affecting optimizations
        return devia * fma(sqrt(Tv(2) / MathConst<Tv>::pi), exp(Tv(-0.5) * square(zscore)), zscore * (Tv(1) + erf(zscore * reciprocal(MathConst<Tv>::sqrt2))));
    }
}
