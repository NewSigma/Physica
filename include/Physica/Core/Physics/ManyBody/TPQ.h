/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Hamilton/HamiltonMatrix.h"

namespace Physica {
    /**
     * NVE ensemble referenced from [1]
     * NVT ensemble referenced from [2]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 108, 240401 (2012); https://doi.org/10.1103/PhysRevLett.108.240401
     * [2] Phys. Rev. Lett. 111, 010401 (2013); https://doi.org/10.1103/PhysRevLett.111.010401
     */
    template<Scalar T>
    class TPQ : public VectorND<T> {
        using This = TPQ<T>;
        using Base = VectorND<T>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;
        // The lowest minus float number, keeping enough digits for later algebra
        constexpr static Tr Smallest = -std::numeric_limits<T>::max() * std::numeric_limits<Tr>::epsilon();

        Tr beta = 0;
        Tr lnZ0 = 0; // Collect constant contribution, improve numerical stability

        Tr traceMu = Tr::nan();
        int numMinCostTerm = 0;
        int numSplit = 0;
    public:
        TPQ() = default;
        TPQ(size_t length);
        TPQ(size_t length, Tr traceMu_, int numMinCostTerm_, int numSplit_);
        TPQ(const This&) = default;
        TPQ(This&&) noexcept = default;
        ~TPQ() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void pre_nvt_step(const Matrix auto& hamiltonH, Tr deltaBeta);
        template<ExecutePolicy P = Sequential>
        void nvt_step(const Matrix auto& hamiltonH, Tr deltaBeta);

        [[nodiscard]] Tr lnPartitionZ() const;
        [[nodiscard]] Tr lnSquaredDot(const VectorND<Tr>& other) const;
        void swap(This& __restrict obj) noexcept;

        template<RNG R = Random<>>
        void random_normal();
        /* Getters */
        [[nodiscard]] Tr getBeta() const noexcept { return beta; }
        [[nodiscard]] Tr getTraceMu() const noexcept;
        [[nodiscard]] int getNumMinCostTerm() const noexcept;
        [[nodiscard]] int getNumSplit() const noexcept;
        [[nodiscard]] const Base& asVector() const noexcept { return *this; }
        [[nodiscard]] Base& asVector() noexcept { return *this; }
        /* Static members */
        template<RNG R = Random<>>
        [[nodiscard]] static This random_normal(size_t len);
        [[nodiscard]] static Tr calcObserveUVT(Tr beta, Tr mu, const MatrixND<T>& lnPartitionNVT, const MatrixND<T>& observeNVT);
    private:
        using Base::random_uniform;
        using Base::random_any;
        [[nodiscard]] bool isPrepared() const noexcept;
        /* Static members */
        constexpr static void checkH(const Matrix auto& hamiltonH) noexcept;
    };

    template<Scalar T>
    TPQ<T>::TPQ(size_t length) : Base(length) {
        assert(length > 0 && "[Error]: Invalid length");
    }

    template<Scalar T>
    TPQ<T>::TPQ(size_t length, Tr traceMu_, int numMinCostTerm_, int numSplit_) : TPQ(length) {
        traceMu = std::move(traceMu_);
        numMinCostTerm = numMinCostTerm_;
        numSplit = numSplit_;
        assert(isPrepared() && "[Error]: Invalid params");
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void TPQ<T>::pre_nvt_step(const Matrix auto& hamiltonH, Tr deltaBeta) {
        checkH(hamiltonH);
        const Tr factor = deltaBeta * Trv(-0.5);
        const auto expr = exp(factor * hamiltonH) * asVector();
        traceMu = expr.calcTraceMu();

        const auto params = expr.template calcParam<P>(traceMu);
        numMinCostTerm = params.first;
        numSplit = params.second;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void TPQ<T>::nvt_step(const Matrix auto& hamiltonH, Tr deltaBeta) {
        checkH(hamiltonH);
        if (deltaBeta.isZero())
            return;

        const Tr factor = deltaBeta * Trv(-0.5);
        const auto expr = exp(factor * hamiltonH) * asVector();
        VectorND<T> dot(Base::getLength());
        if (isPrepared())
            lnZ0 += expr.template assign<true, P>(dot, traceMu, std::make_pair(numMinCostTerm, numSplit));
        else {
            traceMu = expr.calcTraceMu();
            lnZ0 += expr.template assign<true, P>(dot, traceMu);
        }
        Base::swap(dot);
        beta += deltaBeta;
        lnZ0 += traceMu;
    }

    template<Scalar T>
    auto TPQ<T>::lnPartitionZ() const -> Tr {
        if (abs(asVector()).max().isSubNormal())
            return Smallest;
        return Base::lnSquaredNorm() + Tr(2) * lnZ0;
    }

    template<Scalar T>
    auto TPQ<T>::lnSquaredDot(const VectorND<Tr>& other) const -> Tr {
        assert(Base::getLength() == other.getLength() && "[Error]: Dimensions do not match");
        const Tr maxabs = abs(asVector()).max();
        if (maxabs.isSubNormal())
            return Smallest;

        const Tr factor = reciprocal(maxabs);
        const Tr dot = hadamard((asVector() * factor).squaredNorms(), other).sum();
        if (dot.isZero())
            return Smallest;
        return ln(dot) + Tr(2) * (ln(maxabs) + lnZ0);
    }

    template<Scalar T>
    void TPQ<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        beta.swap(obj.beta);
        lnZ0.swap(obj.lnZ0);
        traceMu.swap(obj.traceMu);
        std::swap(numMinCostTerm, obj.numMinCostTerm);
        std::swap(numSplit, obj.numSplit);
    }

    template<Scalar T>
    template<RNG R>
    void TPQ<T>::random_normal() {
        Base::template random_normal<R>();
        beta = 0;
        lnZ0 = 0;
    }

    template<Scalar T>
    template<RNG R>
    TPQ<T> TPQ<T>::random_normal(size_t len) {
        This result(len);
        result.random_normal<R>();
        return result;
    }

    template<Scalar T>
    auto TPQ<T>::getTraceMu() const noexcept -> Tr {
        assert(traceMu.isFinite() && "[Error]: Trace is not ready");
        return traceMu;
    }

    template<Scalar T>
    int TPQ<T>::getNumMinCostTerm() const noexcept {
        assert(isPrepared());
        return numMinCostTerm;
    }

    template<Scalar T>
    int TPQ<T>::getNumSplit() const noexcept {
        assert(isPrepared());
        return numSplit;
    }

    template<Scalar T>
    auto TPQ<T>::calcObserveUVT(Tr beta, Tr mu, const MatrixND<T>& lnPartitionNVT, const MatrixND<T>& observeNVT) -> Tr {
        assert(lnPartitionNVT.getRow() == lnPartitionNVT.getCol());
        assert(lnPartitionNVT.getRow() == observeNVT.getRow());
        assert(lnPartitionNVT.getCol() == observeNVT.getCol());
        const size_t order = observeNVT.getRow();
        auto weights = DenseMatrix<T>(order);
        for (size_t i = 0; i < order; ++i)
            for (size_t j = 0; j < order; ++j)
                weights(i, j) = T(i + j);
        weights = lnPartitionNVT + weights * (beta * mu);
        weights = exp_elem(weights - weights.lnSumExp());
        return hadamard(weights, observeNVT).sum();
    }

    template<Scalar T>
    bool TPQ<T>::isPrepared() const noexcept {
        return traceMu.isFinite() && numMinCostTerm > 0 && numSplit > 0;
    }

    template<Scalar T>
    constexpr void TPQ<T>::checkH(const Matrix auto& hamiltonH) noexcept {
        using M = std::remove_cvref<decltype(hamiltonH)>::type;
        static_assert(std::same_as<T, typename M::ScalarType>, "[Error]: Inconsistent ScalarType");
        static_assert(MatrixOption::isHermiteMatrix<M>(), "[Error]: Do not support non-hermite hamiltionian");
    }
}
