/*
 * Copyright 2020-2024 Weibo He.
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

#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class T> class ComplexScalar;

    namespace Internal {
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<ComplexScalar<AnyScalar1>, AnyScalar2> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = ComplexScalar<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };

        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<AnyScalar1, ComplexScalar<AnyScalar2>> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = ComplexScalar<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };

        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<ComplexScalar<AnyScalar1>, ComplexScalar<AnyScalar2>> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = ComplexScalar<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };
    }

    template<class T>
    class ComplexScalar : public ScalarBase<ComplexScalar<T>> {
        using This = ComplexScalar<T>;
        using Base = ScalarBase<This>;
        using PacketType = typename BestPacket<T, 2>::Type;
    public:
        using typename Base::ScalarType;
        using typename Base::TrivialType;
        constexpr static bool enableSIMD = !std::is_same<T, PacketType>::value;
    private:
        T real;
        T imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(double d);
        ComplexScalar(T real_);
        ComplexScalar(T real_, T imag_);
        ComplexScalar(std::initializer_list<T> list);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept = default;
        /* Operators */
        ComplexScalar& operator=(ComplexScalar c);
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<T>& c) const;
        bool operator!=(const ComplexScalar<T>& c) const { return !(operator==(c)); }
        /* Operations */
        [[nodiscard]] inline PacketType packet() const;
        inline void writePacket(const PacketType packet);
        void swap(ComplexScalar& __restrict c) noexcept;
        /* Getters */
        [[nodiscard]] T squaredNorm() const;
        [[nodiscard]] inline T norm() const;
        [[nodiscard]] T phase() const;
        [[nodiscard]] ComplexScalar unit() const;
        [[nodiscard]] ComplexScalar conjugate() const noexcept { return ComplexScalar(real, -imag); }
        [[nodiscard]] T& getReal() noexcept { return real; }
        [[nodiscard]] const T& getReal() const noexcept { return real; }
        [[nodiscard]] T& getImag() noexcept { return imag; }
        [[nodiscard]] const T& getImag() const noexcept { return imag; }
        [[nodiscard]] bool isZero() const noexcept { return real.isZero() && imag.isZero(); }
        /* Static Members */
        [[nodiscard]] inline static ComplexScalar fromPhase(T phase);
        template<class RandomGenerator>
        [[nodiscard]] static ComplexScalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static ComplexScalar random_normal(RandomGenerator& gen);
        [[nodiscard]] static const H5::DataType& getH5DataType();
    private:
        [[nodiscard]] static H5::DataType* makeH5DataType();
    };

    template<class T>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<T>& c);

    template<class T>
    inline ComplexScalar<T> operator+(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    inline ComplexScalar<T> operator-(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    ComplexScalar<T> operator*(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    ComplexScalar<T> operator/(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const ComplexScalar<ScalarType>& c,const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const ComplexScalar<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const ComplexScalar<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const ComplexScalar<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const Scalar<Option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const Scalar<Option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const Scalar<Option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const Scalar<Option>& s, const ComplexScalar<ScalarType>& c);

    template<class T>
    ComplexScalar<T> operator<<(const ComplexScalar<T>& c, int i) {
        return ComplexScalar<T>(c.getReal() << i, c.getImag() << i);
    }

    template<class T>
    ComplexScalar<T> operator>>(const ComplexScalar<T>& c, int i) {
        return ComplexScalar<T>(c.getReal() >> i, c.getImag() >> i);
    }

    template<class T>
    ComplexScalar<T> operator-(const ComplexScalar<T>& c) { return ComplexScalar(-c.getReal(), -c.getImag()); }

    template<class T>
    void swap(ComplexScalar<T>& __restrict c1, ComplexScalar<T>& __restrict c2) noexcept { c1.swap(c2); }

    template<class T, class ScalarType>
    void operator+=(ComplexScalar<T>& c, const ScalarType& t) { c = c + t; }

    template<class T, class ScalarType>
    void operator-=(ComplexScalar<T>& c, const ScalarType& t) { c = c - t; }

    template<class T, class ScalarType>
    void operator*=(ComplexScalar<T>& c, const ScalarType& t) { c = c * t; }

    template<class T, class ScalarType>
    void operator/=(ComplexScalar<T>& c, const ScalarType& t) { c = c / t; }
}

namespace Physica {
    using namespace Core;

    template<class T>
    class Traits<ComplexScalar<T>> {
        static_assert(!T::isComplex, "[Error]: Double complex mark is not allowed");
        static_assert(!T::isDifferentiable, "[Error]: Differentiable mark should locate in outsite");
    public:
        using ScalarType = ComplexScalar<T>;
        using RealType = T;
        using ComplexType = ComplexScalar<T>;
        using TrivialType = typename T::TrivialType;
        using PlainScalar = ScalarType;
        constexpr static ScalarOption Option = Traits<T>::Option;
        constexpr static bool errorTrack = Traits<T>::errorTrack;
        constexpr static bool isComplex = true;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;
        constexpr static unsigned int Order = 0;
    };
}

namespace std {
    template<class RealType>
    struct numeric_limits<Physica::Core::ComplexScalar<RealType>> : public numeric_limits<RealType> {};
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"
