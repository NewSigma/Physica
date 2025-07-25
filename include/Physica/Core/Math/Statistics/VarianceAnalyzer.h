/*
 * Copyright 2022-2025 Weibo He.
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

namespace Physica {
    /**
     * References:
     * [1] 贾俊平, 何晓群, 金勇进. 统计学.第6版[M]. 中国人民大学出版社, 2015.239-244
     * [2] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016:50-52
     */
    template<Scalar T>
    class VarianceAnalyzer {
        using This = VarianceAnalyzer<T>;
        using DataSet = Array<VectorND<T>>;

        DataSet data;
    public:
        VarianceAnalyzer() = default;
        explicit VarianceAnalyzer(size_t numGroup);
        VarianceAnalyzer(const This&) = default;
        VarianceAnalyzer(This&&) noexcept = default;
        ~VarianceAnalyzer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] VectorND<T>& operator[](size_t i) noexcept { return data[i]; }
        [[nodiscard]] const VectorND<T>& operator[](size_t i) const noexcept { return data[i]; }
        /* Operations */
        [[nodiscard]] T calcRelationCoeff() const noexcept;
        [[nodiscard]] T calcParamF() const;
        [[nodiscard]] T calcParamR() const;
        [[nodiscard]] T calcCorrTime() const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumGroup() const noexcept { return data.getLength(); }
        [[nodiscard]] size_t getTotalNumSample() const noexcept;
    private:
        [[nodiscard]] T calcSSE() const;
        [[nodiscard]] T calcMSE() const;
        [[nodiscard]] T calcSSA() const;
        [[nodiscard]] T calcMSA() const;
        [[nodiscard]] T calcSST() const;
        [[nodiscard]] T calcTotalMean() const;
    };

    template<Scalar T>
    VarianceAnalyzer<T>::VarianceAnalyzer(size_t numGroup) : data(numGroup) {
        assert(numGroup > 1);
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcRelationCoeff() const noexcept {
        T ssa = calcSSA();
        return ssa / (calcSSE() + ssa); // Eq. (10.10) of [1]
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcParamF() const {
        return calcMSA() / calcMSE(); // Eq. (10.9) of [1]
    }
    /**
     * R_X before Eq. (3.20) of [2]
     */
    template<Scalar T>
    T VarianceAnalyzer<T>::calcParamR() const {
        return T(getTotalNumSample() - 1) / T(getNumGroup() - 1) * calcRelationCoeff();
    }
    /**
     * Eq. (3.21) of [2]
     *
     * Note: Param R is required to be saturated
     */
    template<Scalar T>
    T VarianceAnalyzer<T>::calcCorrTime() const {
        return (calcParamR() - 1.0) * 0.5;
    }

    template<Scalar T>
    void VarianceAnalyzer<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        data.swap(obj.data);
    }

    template<Scalar T>
    size_t VarianceAnalyzer<T>::getTotalNumSample() const noexcept {
        size_t count = 0;
        for (const auto& vec : data)
            count += vec.getLength();
        return count;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcSSE() const {
        T result = 0;
        for (const auto& vec : data)
            result += (vec - vec.mean()).squaredNorm();
        return result;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcMSE() const {
        return calcSSE() / T(getTotalNumSample() - getNumGroup());
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcSSA() const {
        const T mean_total = calcTotalMean();
        T result = 0;
        for (const auto& vec : data)
            result += T(vec.getLength()) * square(vec.mean() - mean_total);
        return result;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcMSA() const {
        return calcSSA() / T(getNumGroup() - 1);
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcSST() const {
        return calcSSA() + calcSSE();
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::calcTotalMean() const {
        T result = 0;
        for (const auto& vec : data)
            result += vec.mean() * T(vec.getLength());
        result /= T(getTotalNumSample());
        return result;
    }
}
