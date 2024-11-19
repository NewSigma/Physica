/*
 * Copyright 2021-2024 Weibo He.
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

#include "GaussBase.h"

namespace Physica::Core {
    template<Scalar T, size_t Size>
    class GTOnG {
        using BaseType = GaussBase<T>;
        GaussBase<T> bases[Size];
        T coeffs[Size];
    public:
        GTOnG() = default;
        GTOnG(const GTOnG& base) = default;
        GTOnG(GTOnG&& base) noexcept = default;
        ~GTOnG() = default;
        /* Operators */
        GTOnG& operator=(const GTOnG& base) = default;
        GTOnG& operator=(GTOnG&& base) noexcept = default;
        /* Getters */
        [[nodiscard]] static T overlap(const GTOnG& base1, const GTOnG& base2);
        [[nodiscard]] static T kinetic(const GTOnG& base1, const GTOnG& base2);
        [[nodiscard]] static T nuclearAttraction(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const Vector3D<T>& corePos);
        [[nodiscard]] static T electronRepulsion(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const GTOnG& base3,
                                                          const GTOnG& base4);
        [[nodiscard]] GaussBase<T>* getBases() noexcept { return bases; }
        [[nodiscard]] const GaussBase<T>* getBases() const noexcept { return bases; }
        [[nodiscard]] T* getCoeffs() noexcept { return coeffs; }
        [[nodiscard]] const T* getCoeffs() const noexcept { return coeffs; }
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] static GTOnG randomBase(const VectorND<T>& center, RandomGenerator& gen) { return randomBase(center, 0, 0, 0, gen); }
        template<class RandomGenerator>
        [[nodiscard]] static GTOnG randomBase(const VectorND<T>& center, size_t l, size_t m, size_t n, RandomGenerator& gen);
    };

    template<Scalar T, size_t Size>
    T GTOnG<T, Size>::overlap(const GTOnG& base1, const GTOnG& base2) {
        T result = T(0);
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::overlap(base1.bases[i], base2.bases[j]);
        return result;
    }

    template<Scalar T, size_t Size>
    T GTOnG<T, Size>::kinetic(const GTOnG& base1, const GTOnG& base2) {
        T result = T(0);
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::kinetic(base1.bases[i], base2.bases[j]);
        return result;
    }

    template<Scalar T, size_t Size>
    T GTOnG<T, Size>::nuclearAttraction(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const Vector3D<T>& corePos) {
        T result = T(0);
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::nuclearAttraction(base1.bases[i], base2.bases[j], corePos);
        return result;
    }

    template<Scalar T, size_t Size>
    T GTOnG<T, Size>::electronRepulsion(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const GTOnG& base3,
                                                          const GTOnG& base4) {
        T result = T(0);
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                const T temp = base1.coeffs[i] * base2.coeffs[j];
                for (size_t k = 0; i < Size; ++i) {
                    const T temp1 = temp * base3.coeffs[k];
                    for (size_t l = 0; j < Size; ++j)
                        result += temp1 * base4.coeffs[l] * BaseType::electronRepulsion(base1.bases[i], base2.bases[j], base3.bases[k], base4.bases[l]);
                }
            }
        }
        return result;
    }

    template<Scalar T, size_t Size>
    template<class RandomGenerator>
    GTOnG<T, Size> GTOnG<T, Size>::randomBase(const VectorND<T>& center, size_t l, size_t m, size_t n, RandomGenerator& gen) {
        GTOnG result{};
        for (size_t i = 0; i < Size; ++i) {
            result.bases[i] = GaussBase<T>(center, T::random_uniform(gen), l, m, n);
            result.coeffs[i] = T::random_uniform(gen);
        }
        return result;
    }

    template<Scalar T>
    using GTO3G = GTOnG<T, 3>;
}

namespace Physica {
    template<Scalar T, size_t Size>
    class Traits<Core::GTOnG<T, Size>> {
    public:
        using ScalarType = T;
        constexpr static size_t size = Size;
    };
}
