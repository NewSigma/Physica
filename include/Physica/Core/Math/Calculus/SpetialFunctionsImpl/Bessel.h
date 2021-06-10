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

#include "Physica/Core/Math/Calculus/Chebyshev.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.171
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ0(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T ax = x;
        ax.toAbs();
        T y, ans1, ans2;
        if (ax < T(8)) {
            y = square(ax);
            ans1 = T(57568490574.0) + y * (T(-13362590354.0) + y * (T(651619640.7) + y * (T(-11214424.18) + y * (T(77392.33017) + y * T(-184.9052456)))));
            ans2 = T(57568490411.0) + y * (T(1029532985.0) + y * (T(9494680.718) + y * (T(59272.64853) + y * (T(267.8532712) + y))));
            return ans1 / ans2;
        }
        else {
            const T reciprocal_ax = reciprocal(ax);
            const T z = reciprocal_ax << 3;
            y = square(z);
            ans1 = T(1) + y * (T(-0.1098628627E-2) + y * (T(0.2734510407E-4) + y * (T(-0.2073370639E-5) + y * T(0.2093887211E-6))));
            ans2 = T(-0.1562499995E-1) + y * (T(0.1430488765E-3) + y * (T(-0.6911147651E-5) + y * (T(0.7621095161E-6) - y * T(0.934945152E-7))));
            const T xx = ax - T(0.785398164);
            return sqrt(T(0.636619772) * reciprocal_ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.171
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ1(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T ax = x;
        ax.toAbs();
        T y, ans1, ans2;
        if (ax < T(8)) {
            y = square(ax);
            ans1 = x * (T(72362614232.0) + y * (T(-7895059235.0) + y * (T(242396853.1) + y * (T(-2972611.439) + y * (T(15704.48260) + y * T(-30.16036606))))));
            ans2 = T(144725228442.0) + y * (T(2300535178.0) + y * (T(18583304.74) + y * (T(99447.43394) + y * (T(376.9991397) + y))));
            return ans1 / ans2;
        }
        else {
            const T reciprocal_ax = reciprocal(ax);
            const T z = reciprocal_ax << 3;
            y = square(z);
            ans1 = T(1) + y * (T(0.183105E-2) + y * (T(-0.3516396496E-4) + y * (T(0.2457520174E-5) + y * T(-0.240337019E-6))));
            ans2 = T(0.04687499995) + y * (T(-0.2002690873E-3) + y * (T(0.8449199096E-5) + y * (T(-0.88228987E-6) + y * T(0.105787412E-6))));
            const T xx = ax - T(2.356194491);
            const T result = sqrt(T(0.636619772) * reciprocal_ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
            return x.isNegative() ? -result : result;
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.173
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Integer& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        constexpr int iexp = std::numeric_limits<typename T::TrivialType>::max_exponent / 2;
        constexpr int acc = 160;

        assert(n > 1);
        const T ax = T(x).toAbs();
        const T square_ax = square(ax);
        if (square_ax < T(8 * std::numeric_limits<typename T::TrivialType>::min()))
            return T(0);
        
        const T two_x = T(2) / ax;
        if (ax > T(n)) {
            T bjm = besselJ0(ax); //Possible to optimize: ax is possitive
            T result = besselJ1(ax); //Possible to optimize: ax is possitive
            for (Integer i = 1; i < n; ++i) {
                const T temp = T(i) * two_x * result - bjm;
                bjm = std::move(result);
                result = std::move(temp);
            }
            return result;
        }
        else {
            bool do_sum = false;
            Integer i = n + Integer(sqrt(T(n * acc)));
            i.setByte(0, i.getByte()[0] & (MPUnitMax - 1));
            T bjp(0);
            T bj(1);
            T result(0);
            T sum(0);
            for (; i.isPositive(); --i) {
                T temp = T(i) * two_x * bj - bjp;
                bjp = std::move(bj);
                bj = std::move(temp);
                int k;
                std::frexp(bj.getTrivial(), &k);
                if (k > iexp) {
                    bj = T(std::ldexp(bj.getTrivial(), -iexp));
                    bjp = T(std::ldexp(bjp.getTrivial(), -iexp));
                    result = T(std::ldexp(result.getTrivial(), -iexp));
                    sum = T(std::ldexp(sum.getTrivial(), -iexp));
                }
                if (do_sum) //Optimize: extract loop will avoid branch
                    sum += bj;
                do_sum = !do_sum;
                if (i == n)
                    result = bjp;
            }
            sum = T(2) * sum - bj;
            result /= sum;
            return x.isNegative() && i.isOdd() ? -result : result;
        }
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ(const Integer& n, const Scalar<type, errorTrack>& x) {
        assert(!n.isNegative());
        if (n == 0)
            return besselJ0(x);
        if (n == 1)
            return besselJ1(x);
        return besselJn(n, x);
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.171
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselY0(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T y, ans1, ans2;
        if (x < T(8)) {
            y = square(x);
            ans1 = T(-2957821389.0) + y * (T(7062834065.0) + y * (T(-512359803.6) + y * (T(10879881.29) + y * (T(-86327.92757) + y * T(228.4622733)))));
            ans2 = T(40076544269.0) + y * (T(745249964.8) + y * (T(7189466.438) + y * (T(47447.26470) + y * (T(226.1030244) + y))));
            return (ans1 / ans2) + T(0.636619772) * besselJ0(x) * ln(x);
        }
        else {
            const T reciprocal_x = reciprocal(x);
            const T z = reciprocal_x << 3;
            y = square(z);
            ans1 = T(1) + y * (T(-0.1098628627E-2) + y * (T(0.2734510407E-4) + y * (T(-0.2073370639E-5) + y * T(0.2093887211E-6))));
            ans2 = T(-0.1562499995E-1) + y * (T(0.1430488765E-3) + y * (T(-0.6911147651E-5) + y * (T(0.7621095161E-6) - y * T(0.934945152E-7))));
            const T xx = x - T(0.785398164);
            return sqrt(T(0.636619772) * reciprocal_x) * (sin(xx) * ans1 + z * cos(xx) * ans2);
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.172
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselY1(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T y, ans1, ans2;
        if (x < T(8)) {
            y = square(x);
            ans1 = x * (T(-0.4900604943E13) + y * (T(0.1275274390E13) + y * (T(-0.5153438139E11) + y * (T(0.7349264551E9) + y * (T(-0.4237922726E7) + y * T(0.8511937935E4))))));
            ans2 = T(0.2499580570E14) + y * (T(0.4244419664E12) + y * (T(0.3733650367E10) + y * (T(0.2245904002E8) + y * (T(0.1020426050E6) + y * (T(0.3549632885E3) + y)))));
            return (ans1 / ans2) + T(0.636619772) * (besselJ1(x) * ln(x) - reciprocal(x));
        }
        else {
            const T reciprocal_x = reciprocal(x);
            const T z = reciprocal_x << 3;
            y = square(z);
            ans1 = T(1) + y * (T(-0.183105E-2) + y * (T(-0.3516396496E-4) + y * (T(-0.2457520174E-5) + y * T(0.240337019E-6))));
            ans2 = T(-0.04687499995E-1) + y * (T(-0.2002690873E-3) + y * (T(0.8449199096E-5) + y * (T(-0.88228987E-6) - y * T(0.105787412E-6))));
            const T xx = x - T(2.356194491);
            return sqrt(T(0.636619772) * reciprocal_x) * (sin(xx) * ans1 + z * cos(xx) * ans2);
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.173
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselYn(const Integer& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        assert(n > 1);
        const T two_x = T(2) / x;
        T bym = besselY0(x);
        T result = besselY1(x);
        for (Integer i = 1; i < n; ++i) {
            const T temp = T(i) * two_x * result - bym;
            bym = std::move(result);
            result = std::move(temp);
        }
        return result;
    }

    namespace Internal {
        /**
         * Reference:
         * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
         * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.182
         */
        template<ScalarType type, bool errorTrack>
        void besselChebyshevHelper(
                const Scalar<type, errorTrack>& x
                , Scalar<type, errorTrack>& gamma1
                , Scalar<type, errorTrack>& gamma2
                , Scalar<type, errorTrack>& gamma_plus
                , Scalar<type, errorTrack>& gamma_minus) {
            using T = Scalar<type, errorTrack>;
            assert(abs(x) < T(0.5));
            const static Utils::Array<T> coeff1{-1.142022680371168, 6.5165112670737E-3, 3.087090173086E-4, -3.4706269647E-6, 6.9437664E-9, 3.67795E-11, -1.356E-13};
            const static Utils::Array<T> coeff2{1.843740587300905, -7.68528408447867E-2, 1.2719271366546E-3, -4.9717367042E-6, -3.31261198E-8, 2.423096E-10, -1.702E-13, -1.49E-15};

            const T from(-1);
            const T to(1);
            const T x2 = x * T(2);
            gamma1 = chebyshev_calc_even(from, to, coeff1, x2);
            gamma2 = chebyshev_calc_even(from, to, coeff2, x2);
            gamma_plus = gamma2 - x * gamma1;
            gamma_minus = gamma2 + x * gamma1;
        }
    }
    /**
     * TODO: square of integers is not implemented, using (i * i) instead
     * 
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.180
     */
    template<ScalarType type, bool errorTrack>
    void besselJn_Yn_dJn_dYn(
            const Scalar<type, errorTrack>& n
            , const Scalar<type, errorTrack>& x
            , Scalar<type, errorTrack>& Jn
            , Scalar<type, errorTrack>& Yn
            , Scalar<type, errorTrack>& dJn
            , Scalar<type, errorTrack>& dYn) {
        using T = Scalar<type, errorTrack>;
        constexpr double xmin = 2;
        constexpr double half = 0.5;
        constexpr double pi_trivial = 3.141592653589793;
        constexpr double epsilon_trivial = std::numeric_limits<typename T::TrivialType>::epsilon();
        constexpr double fpmin_trivial = std::numeric_limits<typename T::TrivialType>::min() / epsilon_trivial;
        const T pi = T(pi_trivial);
        const T epsilon = T(epsilon_trivial);
        const T fpmin = T(fpmin_trivial);
        assert(!n.isNegative() && x.isPositive());

        const Integer nl = x < T(xmin)
                             ? Integer(T(n + T(half)))
                             : Integer(T(n - x + T(xmin - half))).isPositive() ? Integer(T(n - x + T(xmin - half))) : Integer(0);
        const T mu = n - T(nl);
        const T square_mu = square(mu);

        const T reciprocal_x = reciprocal(x);
        const T reciprocal_x_2 = reciprocal_x * T(2);
        const T wronskian = reciprocal_x_2 / pi;
        T f_n = n * reciprocal_x;
        if (f_n < fpmin)
            f_n = fpmin;
        bool sign = true;
        /* Lentz method for continued fraction 1 at n */ {
            T b = n * reciprocal_x_2;
            T c(f_n);
            T d(0);
            do {
                b += reciprocal_x_2;
                d = b - d;
                if (abs(d) < fpmin)
                    d = fpmin;
                c = b - reciprocal(c);
                if (abs(c) < fpmin)
                    c = fpmin;
                d = reciprocal(d);
                T delta = c * d;
                f_n *= delta;
                sign = sign ^ d.isNegative(); //If d is negative, sign = -sign
                if (abs(T(delta - T(1))) <= epsilon)
                    break;
            } while(true);
        }
        //get f at \nu = \mu
        const T a_Jv_1 = sign ? fpmin : -fpmin; //a is a constant, 1 is the first iteration
        const T a_dJv_1 = f_n * a_Jv_1;
        T a_Jv_m = a_Jv_1;
        T a_dJv_m = a_dJv_1;
        /* Loop from \nu = n to \nu = \mu */ {
            T factor = n * reciprocal_x;
            for (Integer i = nl - 1; !i.isNegative(); --i) {
                T temp = factor * a_Jv_m + a_dJv_m;
                factor -= reciprocal_x;
                a_dJv_m = factor * temp - a_Jv_m;
                a_Jv_m = std::move(temp);
            }
            if (a_Jv_m.isZero())
                a_Jv_m = epsilon;
        }
        const T f = a_dJv_m / a_Jv_m;

        T Ymu, Ymu_1, Jmu;
        if (x < T(xmin)) { //Temme series method for continued fraction 2
            const T x_2 = x * T(0.5);
            const T pimu = pi * mu;
            const T factor = abs(pimu) < epsilon ? T(1) : pimu / sin(pimu);
            T d = -ln(x_2);
            T e = mu * d;
            const T factor2 = abs(e) < epsilon ? T(1) : sinh(e) / e;
            T gamma1, gamma2, gamma_plus, gamma_minus;
            Internal::besselChebyshevHelper(mu, gamma1, gamma2, gamma_plus, gamma_minus);
            T ff = T(2 / pi_trivial) * factor * (gamma1 * cosh(e) + gamma2 * factor2 * d);
            e = exp(e);
            T p = e / (gamma_plus * pi);
            T q = reciprocal(e * pi * gamma_minus);
            const T pimu_2 = pimu * T(0.5);
            const T factor3 = abs(pimu_2) < epsilon ? T(1) : sin(pimu_2) / pimu_2;
            T r = pi * pimu_2 * square(factor3);
            T c(1);
            d = -square(x_2);
            T sum = ff + r * q;
            T sum1 = p;
            Integer i = 1;
            while (true) {
                const T scalar_i = T(i);
                ff = T(scalar_i * ff + p + q) / (square(scalar_i) - square_mu);
                c *= d / i;
                p /= T(scalar_i - mu);
                q /= T(scalar_i + mu);
                T delta = c * (ff + r * q);
                sum += delta;
                sum1 += T(c * p - scalar_i * delta);
                if (abs(delta) < (T(1) + abs(sum)) * epsilon)
                    break;
                ++i;
            }
            Ymu = -sum;
            Ymu_1 = -sum1 * reciprocal_x_2;
            const T dYmu = mu * reciprocal_x * Ymu - Ymu_1;
            Jmu = wronskian / (dYmu - f * Ymu);
        }
        else { //Lentz method for continued fraction 2
            T a = T(0.25) - square_mu;
            T p = reciprocal_x * T(-0.5);
            T q(1);
            const T br = T(2) * x;
            T bi(2);
            //Start from i = 1
            T factor = a * reciprocal_x / (square(p) + square(q));
            T cr = br + q * factor;
            T ci = bi + p * factor;
            T den = square(br) + square(bi);
            T dr = br / den;
            T di = -bi / den;
            T dlr = cr * dr - ci * di;
            T dli = cr * di + ci * dr;
            T temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            Integer i = 1;
            while (true) {
                a += T(i << 1U);
                bi += T(2);
                dr = a * dr + br;
                di = a * di + bi;
                if ((abs(dr) + abs(di)) < fpmin)
                    dr = fpmin;
                factor = a / (square(cr) + square(ci));
                cr = br + cr * factor;
                ci = bi - ci * factor;
                if ((abs(cr) + abs(ci)) < fpmin)
                    cr = fpmin;
                den = square(dr) + square(di);
                dr /= den;
                di /= -den;
                dlr = cr * dr - ci * di;
                dli = cr * di + ci * dr;
                temp = p * dlr - q * dli;
                q = p * dli + q * dlr;
                p = temp;
                if ((abs(T(dlr - T(1))) + abs(dli)) < epsilon)
                    break;
                ++i;
            }
            const T gamma = (p - f) / q;
            Jmu = sqrt(wronskian / ((p - f) * gamma + q));
            Jmu = a_Jv_m > 0 ? Jmu : -Jmu;
            Ymu = Jmu * gamma;
            const T dYmu = Ymu * (p + q / gamma);
            Ymu_1 = mu * reciprocal_x * Ymu - dYmu;
        }
        const T factor = Jmu / a_Jv_m;
        Jn = a_Jv_1 * factor;
        dJn = a_dJv_1 * factor;
        //Loop from v = \mu to v = n
        T Yv = Ymu, Yv_1 = Ymu_1;
        for (Integer i = 1; i <= nl; ++i) {
            T temp = (mu + T(i)) * reciprocal_x_2 * Yv_1 - Yv;
            Yv = Yv_1;
            Yv_1 = std::move(temp);
        }
        Yn = Yv;
        dYn = n * reciprocal_x * Yv - Yv_1;
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return Jn;
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besseldJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return dJn;
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return Yn;
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besseldYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return dYn;
    }
}
