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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
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

        Tr beta;
        Tr traceMu;
        int numMinCostTerm;
        int numSplit;
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
        template<class ModelType, class Executor = SeqExecutor>
        void pre_nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, Tr deltaBeta);
        template<class ModelType, class Executor = SeqExecutor>
        void nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, Tr deltaBeta);

        [[nodiscard]] inline Tr calcPartitionXi() const;
        [[nodiscard]] inline Tr lnPartitionXi() const;
        [[nodiscard]] Tr lnSquaredDot(const VectorND<Tr>& other) const;
        void swap(This& __restrict obj) noexcept;

        template<RNG R>
        inline void random_normal(Tr norm = 0);
        /* Getters */
        [[nodiscard]] Tr getBeta() const noexcept { return beta; }
        [[nodiscard]] Tr getTraceMu() const noexcept { return traceMu; }
        [[nodiscard]] int getNumMinCostTerm() const noexcept { return numMinCostTerm; }
        [[nodiscard]] int getNumSplit() const noexcept { return numSplit; }
        [[nodiscard]] const Base& asVector() const noexcept { return *this; }
        [[nodiscard]] Base& asVector() noexcept { return *this; }
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t len, Tr norm);
    private:
        using Base::random_uniform;
        using Base::random_any;
        [[nodiscard]] bool isPrepared() const;
    };

    template<Scalar T>
    TPQ<T>::TPQ(size_t length)
            : Base(length)
            , beta(0)
            , traceMu(std::numeric_limits<Tr>::max())
            , numMinCostTerm(0)
            , numSplit(0) {
        assert(length > 0 && "[Error]: Invalid length");
    }

    template<Scalar T>
    TPQ<T>::TPQ(size_t length, Tr traceMu_, int numMinCostTerm_, int numSplit_)
            : TPQ(length) {
        traceMu = std::move(traceMu_);
        numMinCostTerm = numMinCostTerm_;
        numSplit = numSplit_;
        assert(isPrepared() && "[Error]: Invalid params");
    }

    template<Scalar T>
    template<class ModelType, class Executor>
    void TPQ<T>::pre_nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, Tr deltaBeta) {
        const Tr factor = deltaBeta * Trv(-0.5);
        const auto& hamiltonH = hamiltonH_.getDerived();
        const auto expr = exp(factor * hamiltonH) * asVector();
        traceMu = expr.calcTraceMu();
        const auto params = expr.template calcParam<Executor>(traceMu);
        numMinCostTerm = params.first;
        numSplit = params.second;
    }

    template<Scalar T>
    template<class ModelType, class Executor>
    void TPQ<T>::nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, Tr deltaBeta) {
        if (deltaBeta.isZero())
            return;
        using BufferType = VectorND<T>;
        const Tr factor = deltaBeta * Trv(-0.5);
        const auto& hamiltonH = hamiltonH_.getDerived();
        const auto expr = exp(factor * hamiltonH) * asVector();
        BufferType dot(Base::getLength());
        if (isPrepared())
            expr.template assign<BufferType, Executor>(dot, traceMu, std::make_pair(numMinCostTerm, numSplit));
        else
            expr.template assign<BufferType, Executor>(dot);
        Base::swap(dot);
        beta += deltaBeta;
    }

    template<Scalar T>
    auto TPQ<T>::calcPartitionXi() const -> Tr {
        return Base::squaredNorm();
    }

    template<Scalar T>
    auto TPQ<T>::lnPartitionXi() const -> Tr {
        const bool isUnderflow = abs(asVector()).max() < Tr(std::numeric_limits<Tr>::min());
        if (isUnderflow)
            return Tr(-std::numeric_limits<T>::max());
        return Base::lnSquaredNorm();
    }

    template<Scalar T>
    auto TPQ<T>::lnSquaredDot(const VectorND<Tr>& other) const -> Tr {
        assert(Base::getLength() == other.getLength() && "[Error]: Dimensions do not match");
        const Tr maxabs = abs(asVector()).max();
        const bool isUnderflow = maxabs < Tr(std::numeric_limits<Tr>::min());
        if (isUnderflow)
            return Tr(-std::numeric_limits<T>::max());
        const Tr factor = reciprocal(maxabs);
        const auto expr1 = asVector() * factor;
        const auto expr2 = expr1.squaredNorms();
        const Tr dot = hadamard(expr2, other).sum();
        if (dot.isZero())
            return Tr(-std::numeric_limits<T>::max());
        return ln(dot) + Tr(2) * ln(maxabs);
    }

    template<Scalar T>
    void TPQ<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        beta.swap(obj.beta);
        traceMu.swap(obj.traceMu);
        std::swap(numMinCostTerm, obj.numMinCostTerm);
        std::swap(numSplit, obj.numSplit);
    }

    template<Scalar T>
    template<RNG R>
    inline void TPQ<T>::random_normal(Tr norm) {
        const bool useDefault = norm.isZero();
        if (useDefault) // Default norm is as small as possible while keep all effective digits
            norm = Tr(std::numeric_limits<Tr>::min() / std::numeric_limits<Tr>::epsilon());
        Base::template random_normal<R>();
        asVector() *= norm;
    }

    template<Scalar T>
    template<RNG R>
    TPQ<T> TPQ<T>::random_normal(size_t len, Tr norm) {
        This result(len);
        result.random_normal<R>(norm);
        return result;
    }

    template<Scalar T>
    bool TPQ<T>::isPrepared() const {
        const bool muReady = traceMu != Tr(std::numeric_limits<Tr>::max());
        return muReady && numMinCostTerm > 0 && numSplit > 0;
    }
}
