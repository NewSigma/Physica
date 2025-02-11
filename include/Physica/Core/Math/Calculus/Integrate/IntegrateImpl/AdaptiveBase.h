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

#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica {
    template<Scalar T, bool TakeLn>
    class AdaptiveBase {
        using This = AdaptiveBase<T, TakeLn>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;
    protected:
        VectorND<Trv> from;
        VectorND<Trv> to;
        VectorND<T> means;
        VectorND<T> vars;
        VectorND<Trv> loss;
    private:
        int numRefine;
        int numSample;
    public:
        AdaptiveBase(VectorND<Trv> from_, VectorND<Trv> to_, int numRefine_, int numSample_);
        AdaptiveBase(const This&) = default;
        AdaptiveBase(This&&) noexcept = default;
        ~AdaptiveBase() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calcMean(int from = 0) const requires(!TakeLn);
        [[nodiscard]] T calcDevia(int from = 0) const requires(!TakeLn);
        [[nodiscard]] T calcVar(int from = 0) const requires(!TakeLn);

        [[nodiscard]] T calcLnMean(int from = 0) const requires(TakeLn);
        [[nodiscard]] T calcLnDevia(int from = 0) const requires(TakeLn);
        [[nodiscard]] T calcLnVar(int from = 0) const requires(TakeLn);

        [[nodiscard]] T calcSquaredChi(int from = 0) const;

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return from.getLength(); }
        [[nodiscard]] int getNumRefine() const noexcept { return numRefine; }
        [[nodiscard]] int getNumSample() const noexcept { return numSample; }
        [[nodiscard]] const auto& getMeans() const noexcept { return means; }
        [[nodiscard]] const auto& getVars() const noexcept { return vars; }
        [[nodiscard]] const auto& getLoss() const noexcept { return loss; }
        /* Setters */
        void setNumRefine(int numRefine_);
    };

    template<Scalar T, bool TakeLn>
    AdaptiveBase<T, TakeLn>::AdaptiveBase(VectorND<Trv> from_, VectorND<Trv> to_, int numRefine_, int numSample_)
            : from(std::move(from_))
            , to(std::move(to_))
            , numSample(numSample_) {
        assert(from.getLength() == to.getLength() && "[Error]: Inconsistent dim");
        assert(numSample > 0);

        setNumRefine(numRefine_);
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcMean(int from) const requires(!TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return reciprocal(vars.tail(from)) * means.tail(from) / reciprocal(vars.tail(from)).sum();
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcVar(int from) const requires(!TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return reciprocal(reciprocal(vars.tail(from)).sum());
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcDevia(int from) const requires(!TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return sqrt(calcVar(from));
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcLnMean(int from) const requires(TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return (means.tail(from) - vars.tail(from)).lnSumExp() + calcLnVar();
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcLnDevia(int from) const requires(TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return Trv(0.5) * calcLnVar(from);
    }

    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcLnVar(int from) const requires(TakeLn) {
        assert(0 <= from && from < getNumRefine());
        return -(-vars.tail(from)).lnSumExp();
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:414-416
     */
    template<Scalar T, bool TakeLn>
    T AdaptiveBase<T, TakeLn>::calcSquaredChi(int from) const {
        assert(0 <= from && from < getNumRefine());
        if (numRefine == 1)
            return 1;
        if constexpr (TakeLn) {
            const T lnMean = calcLnMean(from);
            VectorND<T> buffer = exp(means.tail(from) - lnMean);
            buffer = ln(abs(buffer - Trv(1))) + lnMean;
            buffer = Trv(2) * buffer - vars.tail(from);
            return exp(buffer).sum() / Trv(numRefine - from - 1);
        }
        else {
            const Trv factor = Trv(numSample) / Trv(numRefine - from - 1);
            const T mean1 = calcMean(from);
            return factor * divide((means.tail(from) - mean1).squaredNorms(), vars.tail(from)).sum(); // Normalize, refer to [1]
        }
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, bool TakeLn>
    const H5Group AdaptiveBase<T, TakeLn>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        group.readAttr("NumRefine", numRefine);
        group.readAttr("NumSample", numSample);

        from.read(group, "From");
        to.read(group, "To");
        means.read(group, "Means");
        vars.read(group, "Vars");
        loss.read(group, "Loss");
        return group;
    }

    template<Scalar T, bool TakeLn>
    H5Group AdaptiveBase<T, TakeLn>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        group.writeAttr("NumRefine", numRefine);
        group.writeAttr("NumSample", numSample);

        from.write(group, "From");
        to.write(group, "To");
        means.write(group, "Means");
        vars.write(group, "Vars");
        loss.write(group, "Loss");
        return group;
    }
#endif

    template<Scalar T, bool TakeLn>
    void AdaptiveBase<T, TakeLn>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        from.swap(obj.from);
        to.swap(obj.to);
        std::swap(numRefine, obj.numRefine);
        std::swap(numSample, obj.numSample);

        means.swap(obj.means);
        vars.swap(obj.vars);
        loss.swap(obj.loss);
    }

    template<Scalar T, bool TakeLn>
    void AdaptiveBase<T, TakeLn>::setNumRefine(int numRefine_) {
        assert(numRefine_ > 0);
        numRefine = numRefine_;
        means.resize(numRefine);
        vars.resize(numRefine);
        loss.resize(numRefine);
        means = Trv(0);
        vars = Trv(0);
    }
}
