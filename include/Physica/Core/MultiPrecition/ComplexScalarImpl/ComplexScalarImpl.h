/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXSCALARIMPL_H
#define PHYSICA_COMPLEXSCALARIMPL_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 */
namespace Physica::Core {
    //!Optimize: maybe use && to avoid unnecessary move.
    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>::ComplexScalar(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2)
            : real(std::move(s1)), imagine(std::move(s2)) {}

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>::ComplexScalar(ComplexScalar&& c) noexcept
            : real(std::move(c.real)), imagine(std::move(c.imagine)) { Q_UNUSED(type) Q_UNUSED(errorTrack) }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>& ComplexScalar<type, errorTrack>::operator=(const ComplexScalar& c) {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        if(this == &c)
            return *this;
        real = c.real;
        imagine = c.imagine;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>& ComplexScalar<type, errorTrack>::operator=(ComplexScalar&& c) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        real = std::move(c.real);
        imagine = std::move(c.imagine);
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    void ComplexScalar<type, errorTrack>::swap(ComplexScalar& c) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        Physica::Core::swap(real, c.real);
        Physica::Core::swap(imagine, c.imagine);
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getZero() {
        return ComplexScalar(Scalar<type, errorTrack>::getZero(), Scalar<type, errorTrack>::getZero());
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getOne() {
        return ComplexScalar(Scalar<type, errorTrack>::getOne(), Scalar<type, errorTrack>::getZero());
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getRandom() {
        return ComplexScalar(randomScalar<type, errorTrack>(), randomScalar<type, errorTrack>());
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> norm(const ComplexScalar<type, errorTrack>& c) {
        return sqrt(square(c.getReal()) + square(c.getImagine()));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arg(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine  = c.getImagine();
        auto result = arctan(imagine / real);
        //arctan is defined on [-Pi / 2, Pi / 2], judging the quadrant is necessary.
        if(real.isPositive()) {
            if(imagine.isNegative())
                result += MathConst::getInstance().getPI() << 1;
        }
        else
            result += MathConst::getInstance().getPI();
        return result;
    }

    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<type, errorTrack>& c) {
        const auto& imagine = c.getImagine();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : "+" )<< 'i' << fabs(double(imagine)) << std::setprecision(6);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator+(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        return ComplexScalar<type, errorTrack>(c1.getReal() + c2.getReal(), c1.getImagine() + c2.getImagine());
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator-(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        return ComplexScalar<type, errorTrack>(c1.getReal() - c2.getReal(), c1.getImagine() - c2.getImagine());
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator*(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImagine();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImagine();
        /*
         * Optimize:
         * Use (a + ib)(c + id) = (ac - bd) + i((a + b)(c + d) - ac - bd)
         * instead of (a + ib)(c + id) = (ac - bd) + i(ad + bc) to avoid multiply.
         * But it is unclear if this method is useful to every machine.
         * May be add checks and use Config.h to determine which method to use.
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        return ComplexScalar<type, errorTrack>(ac - bd
                , (real_1 + real_2) * (imagine_1 + imagine_2) - ac - bd);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator/(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImagine();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImagine();
        /*
         * Optimize: Using the same method with operator*().
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        const auto divisor = square(real_2) + square(imagine_2);
        return ComplexScalar<type, errorTrack>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() + s, c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() - s, c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() * s, c.getImagine() * s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() / s, c.getImagine() / s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() + s, c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(s - c.getReal(), c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() * s, c.getImagine() * s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImagine();
        const auto divisor = s * reciprocal(square(real) + square(imagine));
        return ComplexScalar<type, errorTrack1 | errorTrack2>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"

#endif