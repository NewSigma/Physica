/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core::Physics {
    template<class ScalarType, size_t Size>
    class GTOnG;

    namespace Internal {
        template<class T> class Traits;

        template<class T, size_t Size>
        class Traits<GTOnG<T, Size>> {
        public:
            using ScalarType = T;
            constexpr static size_t size = Size;
        };
    }

    template<class ScalarType, size_t Size>
    class GTOnG {
        using BaseType = GaussBase<ScalarType>;
        GaussBase<ScalarType> bases[Size];
        ScalarType coeffs[Size];
    public:
        GTOnG() = default;
        GTOnG(const GTOnG& base) = default;
        GTOnG(GTOnG&& base) noexcept = default;
        ~GTOnG() = default;
        /* Operators */
        GTOnG& operator=(const GTOnG& base) = default;
        GTOnG& operator=(GTOnG&& base) noexcept = default;
        /* Getters */
        [[nodiscard]] static ScalarType overlap(const GTOnG& base1, const GTOnG& base2);
        [[nodiscard]] static ScalarType kinetic(const GTOnG& base1, const GTOnG& base2);
        [[nodiscard]] static ScalarType nuclearAttraction(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const Vector<ScalarType, 3>& corePos);
        [[nodiscard]] static ScalarType electronRepulsion(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const GTOnG& base3,
                                                          const GTOnG& base4);
        [[nodiscard]] GaussBase<ScalarType>* getBases() noexcept { return bases; }
        [[nodiscard]] const GaussBase<ScalarType>* getBases() const noexcept { return bases; }
        [[nodiscard]] ScalarType* getCoeffs() noexcept { return coeffs; }
        [[nodiscard]] const ScalarType* getCoeffs() const noexcept { return coeffs; }
        /* Static Members */
        [[nodiscard]] static GTOnG randomBase(const Vector<ScalarType>& center) { return randomBase(center, 0, 0, 0); }
        [[nodiscard]] static GTOnG randomBase(const Vector<ScalarType>& center, size_t l, size_t m, size_t n);
    };

    template<class ScalarType, size_t Size>
    ScalarType GTOnG<ScalarType, Size>::overlap(const GTOnG& base1, const GTOnG& base2) {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::overlap(base1.bases[i], base2.bases[j]);
        return result;
    }

    template<class ScalarType, size_t Size>
    ScalarType GTOnG<ScalarType, Size>::kinetic(const GTOnG& base1, const GTOnG& base2) {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::kinetic(base1.bases[i], base2.bases[j]);
        return result;
    }

    template<class ScalarType, size_t Size>
    ScalarType GTOnG<ScalarType, Size>::nuclearAttraction(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const Vector<ScalarType, 3>& corePos) {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < Size; ++i)
            for (size_t j = 0; j < Size; ++j)
                result += base1.coeffs[i] * base2.coeffs[j] * BaseType::nuclearAttraction(base1.bases[i], base2.bases[j], corePos);
        return result;
    }

    template<class ScalarType, size_t Size>
    ScalarType GTOnG<ScalarType, Size>::electronRepulsion(const GTOnG& base1,
                                                          const GTOnG& base2,
                                                          const GTOnG& base3,
                                                          const GTOnG& base4) {
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                const ScalarType temp = base1.coeffs[i] * base2.coeffs[j];
                for (size_t k = 0; i < Size; ++i) {
                    const ScalarType temp1 = temp * base3.coeffs[k];
                    for (size_t l = 0; j < Size; ++j)
                        result += temp1 * base4.coeffs[l] * BaseType::electronRepulsion(base1.bases[i], base2.bases[j], base3.bases[k], base4.bases[l]);
                }
            }
        }
        return result;
    }

    template<class ScalarType, size_t Size>
    GTOnG<ScalarType, Size> GTOnG<ScalarType, Size>::randomBase(const Vector<ScalarType>& center, size_t l, size_t m, size_t n) {
        GTOnG result{};
        for (size_t i = 0; i < Size; ++i) {
            result.bases[i] = GaussBase<ScalarType>(center, randomScalar<ScalarType>(), l, m, n);
            result.coeffs[i] = randomScalar<ScalarType>();
        }
        return result;
    }

    template<class ScalarType>
    using GTO3G = GTOnG<ScalarType, 3>;
}
