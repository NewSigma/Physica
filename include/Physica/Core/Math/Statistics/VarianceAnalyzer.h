/*
 * Copyright 2022-2026 Weibo He.
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

#include <bit>
#include <optional>
#include <ranges>
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
        [[nodiscard]] Vector2D<T> statistics() const;
        [[nodiscard]] T calcRelationCoeff() const noexcept;
        [[nodiscard]] T calcParamF() const;
        [[nodiscard]] T calcParamR() const;
        [[nodiscard]] T calcCorrTime() const;

        [[nodiscard]] constexpr auto begin(this auto&& self) noexcept { return self.data.begin(); }
        [[nodiscard]] constexpr auto end(this auto&& self) noexcept { return self.data.end(); }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t size() const noexcept { return getNumGroup(); }
        [[nodiscard]] size_t getNumGroup() const noexcept { return data.getLength(); }
        [[nodiscard]] size_t getTotalNumSample() const noexcept;
        /* Static members */
        [[nodiscard]] static auto blocking(const VectorND<T>& samples) -> std::optional<This>;
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
    auto VarianceAnalyzer<T>::statistics() const -> Vector2D<T> {
        return {calcTotalMean(), sqrt(calcMSA() / T(getTotalNumSample()))};
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
     * Add 0.5 to the expression in Eq. (3.21) of [2] to align with the convention of \class Correlation
     *
     * Note: Param R is required to be saturated
     */
    template<Scalar T>
    T VarianceAnalyzer<T>::calcCorrTime() const {
        return calcParamR() * T(0.5);
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
        assert(!data.empty());
        T result = 0;
        for (const auto& vec : data)
            result += vec.mean() * T(vec.getLength());
        result /= T(getTotalNumSample());
        return result;
    }
    /**
     * \param samples a autocorrelated sequence
     */
    template<Scalar T>
    auto VarianceAnalyzer<T>::blocking(const VectorND<T>& samples) -> std::optional<This> {
        constexpr size_t MinBlockSize = 2;
        struct Config {
            size_t blockSize;
            size_t numGroup;
            T error;
        };

        const size_t numSample = samples.getLength();
        assert(numSample > 1 && "[Error]: At least two samples are required");
        const auto configs = Array<Config>::generate([&](size_t i) {
            const size_t m = MinBlockSize << i;
            const size_t numGroup = numSample / m;
            VectorND<T> means(numGroup);
            for (const size_t g : std::views::iota(0UZ, numGroup))
                means[g] = samples.segment(g * m, (g + 1) * m).mean();
            return Config{m, numGroup, means.deviation(means.mean()) / sqrt(T(numGroup))};
        }, std::bit_width(numSample / (2 * MinBlockSize)));

        const T z95 = T(1.96);
        const T certify = T(0.3);
        for (size_t i = 0; i < configs.getLength(); ++i) {
            const auto& proposal = configs[i];
            bool obvious = true;
            bool verified = false;
            for (size_t j = i + 1; j < configs.getLength(); ++j) {
                const auto& config = configs[j];
                if (proposal.error.isSubNormal()) [[unlikely]]
                    obvious = config.error.isSubNormal();
                else {
                    const T fine = T(proposal.numGroup);
                    const T coarse = T(config.numGroup);
                    const T band = z95 * MathConst<T>::sqrt2 * sqrt((fine - coarse) / ((fine - T(1)) * (coarse - T(1))));
                    bool uncertain = band > certify;
                    if (uncertain)
                        continue;
                    obvious = abs(square(config.error / proposal.error) - T(1)) <= band;
                }

                verified = true;
                if (!obvious)
                    break;
            }
            if (!verified || !obvious)
                continue;

            const size_t numGroup = proposal.numGroup;
            This result(numGroup);
            for (size_t g = 0; g < numGroup; ++g)
                result[g] = samples.segment(g * proposal.blockSize, (g + 1) * proposal.blockSize);
            return result;
        }
        return std::nullopt;
    }
}
