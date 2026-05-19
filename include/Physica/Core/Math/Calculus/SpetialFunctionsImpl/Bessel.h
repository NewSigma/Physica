/*
 * Copyright 2021-2026 Weibo He.
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

namespace Physica {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:171
     */
    template<FloatPrec Prec>
    Real<Prec> besselJ0(const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T ax = abs(x);
        if (ax < T(8)) {
            T y = square(ax);
            T ans1 = fma(y, fma(y, fma(y, fma(y, fma(y, T(-184.9052456), T(77392.33017)), T(-11214424.18)), T(651619640.7)), T(-13362590354.0)), T(57568490574.0));
            T ans2 = fma(y, fma(y, fma(y, fma(y, y + T(267.8532712), T(59272.64853)), T(9494680.718)), T(1029532985.0)), T(57568490411.0));
            return ans1 / ans2;
        }

        const T reciprocal_ax = reciprocal(ax);
        const T z = reciprocal_ax << 3;
        T y = square(z);
        T ans1 = fma(y, fma(y, fma(y, fma(y, T(0.2093887211E-6), T(-0.2073370639E-5)), T(0.2734510407E-4)), T(-0.1098628627E-2)), T(1));
        T ans2 = fma(y, fma(y, fma(y, fma(y, -T(0.934945152E-7), T(0.7621095161E-6)), T(-0.6911147651E-5)), T(0.1430488765E-3)), T(-0.1562499995E-1));
        const T xx = ax - T(0.785398164);
        return sqrt(T(0.636619772) * reciprocal_ax) * fma(cos(xx), ans1, -z * sin(xx) * ans2);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:171
     */
    template<FloatPrec Prec>
    Real<Prec> besselJ1(const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T ax = abs(x);
        T y, ans1, ans2;
        if (ax < T(8)) {
            y = square(ax);
            ans1 = x * fma(y, fma(y, fma(y, fma(y, fma(y, T(-30.16036606), T(15704.48260)), T(-2972611.439)), T(242396853.1)), T(-7895059235.0)), T(72362614232.0));
            ans2 = fma(y, fma(y, fma(y, fma(y, y + T(376.9991397), T(99447.43394)), T(18583304.74)), T(2300535178.0)), T(144725228442.0));
            return ans1 / ans2;
        }

        const T reciprocal_ax = reciprocal(ax);
        const T z = reciprocal_ax << 3;
        y = square(z);
        ans1 = fma(y, fma(y, fma(y, fma(y, T(-0.240337019E-6), T(0.2457520174E-5)), T(-0.3516396496E-4)), T(0.183105E-2)), T(1));
        ans2 = fma(y, fma(y, fma(y, fma(y, T(0.105787412E-6), T(-0.88228987E-6)), T(0.8449199096E-5)), T(-0.2002690873E-3)), T(0.04687499995));
        const T xx = ax - T(2.356194491);
        const T result = sqrt(T(0.636619772) * reciprocal_ax) * fma(cos(xx), ans1, -z * sin(xx) * ans2);
        return x.isNegative() ? -result : result;
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:173
     */
    template<FloatPrec Prec>
    Real<Prec> besselJn(const Integer& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        constexpr int iexp = std::numeric_limits<typename T::MachineType>::max_exponent / 2;
        constexpr int acc = 160;

        assert(n > 1);
        const T ax = abs(x);
        const T square_ax = square(ax);
        if (square_ax < T(8 * std::numeric_limits<typename T::MachineType>::min()))
            return T(0);

        const T two_x = T(2) / ax;
        if (ax > T(n)) {
            T bjm = besselJ0(ax); // Possible to optimize: ax is possitive
            T result = besselJ1(ax); // Possible to optimize: ax is possitive
            for (Integer i = 1; i < n; ++i) {
                const T temp = fma(T(i) * two_x, result, -bjm);
                bjm = std::move(result);
                result = std::move(temp);
            }
            return result;
        }

        bool do_sum = false;
        Integer i = n + Integer(sqrt(T(n * acc)));
        i.setByte(0, i.getByte()[0] & (MPUnitMax - 1));
        T bjp(0);
        T bj(1);
        T result(0);
        T sum(0);
        for (; i.isPositive(); --i) {
            T temp = fma(T(i) * two_x, bj, -bjp);
            bjp = std::move(bj);
            bj = std::move(temp);
            int k;
            std::frexp(bj.toMachine(), &k);
            if (k > iexp) {
                bj = T(std::ldexp(bj.toMachine(), -iexp));
                bjp = T(std::ldexp(bjp.toMachine(), -iexp));
                result = T(std::ldexp(result.toMachine(), -iexp));
                sum = T(std::ldexp(sum.toMachine(), -iexp));
            }
            if (do_sum) // Optimize: extract loop will avoid branch
                sum += bj;
            do_sum = !do_sum;
            if (i == n)
                result = bjp;
        }
        sum = fma(T(2), sum, -bj);
        result /= sum;
        return x.isNegative() && i.isOdd() ? -result : result;
    }

    template<FloatPrec Prec>
    Real<Prec> besselJ(const Integer& n, const Real<Prec>& x) noexcept {
        assert(!n.isNegative());
        if (n == 0)
            return besselJ0(x);
        if (n == 1)
            return besselJ1(x);
        return besselJn(n, x);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:171
     */
    template<FloatPrec Prec>
    Real<Prec> besselY0(const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        if (x < T(8)) {
            T y = square(x);
            T ans1 = fma(y, fma(y, fma(y, fma(y, fma(y, T(228.4622733), T(-86327.92757)), T(10879881.29)), T(-512359803.6)), T(7062834065.0)), T(-2957821389.0));
            T ans2 = fma(y, fma(y, fma(y, fma(y, y + T(226.1030244), T(47447.26470)), T(7189466.438)), T(745249964.8)), T(40076544269.0));
            return fma(T(0.636619772) * besselJ0(x), ln(x), ans1 / ans2);
        }
        const T reciprocal_x = reciprocal(x);
        const T z = reciprocal_x << 3;
        T y = square(z);
        T ans1 = fma(y, fma(y, fma(y, fma(y, T(0.2093887211E-6), T(-0.2073370639E-5)), T(0.2734510407E-4)), T(-0.1098628627E-2)), T(1));
        T ans2 = fma(y, fma(y, fma(y, fma(y, -T(0.934945152E-7), T(0.7621095161E-6)), T(-0.6911147651E-5)), T(0.1430488765E-3)), T(-0.1562499995E-1));
        const T xx = x - T(0.785398164);
        return sqrt(T(0.636619772) * reciprocal_x) * fma(sin(xx), ans1, z * cos(xx) * ans2);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:172
     */
    template<FloatPrec Prec>
    Real<Prec> besselY1(const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        if (x < T(8)) {
            T y = square(x);
            T ans1 = x * fma(y, fma(y, fma(y, fma(y, fma(y, T(0.8511937935E4), T(-0.4237922726E7)), T(0.7349264551E9)), T(-0.5153438139E11)), T(0.1275274390E13)), T(-0.4900604943E13));
            T ans2 = fma(y, fma(y, fma(y, fma(y, fma(y, y + T(0.3549632885E3), T(0.1020426050E6)), T(0.2245904002E8)), T(0.3733650367E10)), T(0.4244419664E12)), T(0.2499580570E14));
            return fma(T(0.636619772), besselJ1(x) * ln(x) - reciprocal(x), ans1 / ans2);
        }
        const T reciprocal_x = reciprocal(x);
        const T z = reciprocal_x << 3;
        T y = square(z);
        T ans1 = fma(y, fma(y, fma(y, fma(y, T(0.240337019E-6), T(-0.2457520174E-5)), T(-0.3516396496E-4)), T(-0.183105E-2)), T(1));
        T ans2 = fma(y, fma(y, fma(y, fma(y, -T(0.105787412E-6), T(-0.88228987E-6)), T(0.8449199096E-5)), T(-0.2002690873E-3)), T(-0.04687499995E-1));
        const T xx = x - T(2.356194491);
        return sqrt(T(0.636619772) * reciprocal_x) * fma(sin(xx), ans1, z * cos(xx) * ans2);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:173
     */
    template<FloatPrec Prec>
    Real<Prec> besselYn(const Integer& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        assert(n > 1);
        const T two_x = T(2) / x;
        T bym = besselY0(x);
        T result = besselY1(x);
        for (Integer i = 1; i < n; ++i) {
            const T temp = fma(T(i) * two_x, result, -bym);
            bym = std::move(result);
            result = std::move(temp);
        }
        return result;
    }

    namespace Internal {
        /**
         * Reference:
         * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:182
         */
        template<FloatPrec Prec>
        void besselChebyshevHelper(const Real<Prec>& x, Real<Prec>& gamma1, Real<Prec>& gamma2, Real<Prec>& gamma_plus, Real<Prec>& gamma_minus) noexcept {
            using T = Real<Prec>;
            assert(abs(x) <= T(0.5));
            const static DenseVector<T, 7> coeff1{-1.142022680371168, 6.5165112670737E-3, 3.087090173086E-4, -3.4706269647E-6, 6.9437664E-9, 3.67795E-11, -1.356E-13};
            const static DenseVector<T, 8> coeff2{1.843740587300905, -7.68528408447867E-2, 1.2719271366546E-3, -4.9717367042E-6, -3.31261198E-8, 2.423096E-10, -1.702E-13, -1.49E-15};

            const T from(-1);
            const T to(1);
            const T x2 = x * T(2);
            gamma1 = chebyshev_calc_even(from, to, coeff1, x2);
            gamma2 = chebyshev_calc_even(from, to, coeff2, x2);
            gamma_plus = fma(-x, gamma1, gamma2);
            gamma_minus = fma(x, gamma1, gamma2);
        }
    }
    /**
     * TODO: square of integers is not implemented, using (i * i) instead
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:180
     */
    template<FloatPrec Prec>
    void besselJn_Yn_dJn_dYn(const Real<Prec>& n, const Real<Prec>& x, Real<Prec>& Jn, Real<Prec>& Yn, Real<Prec>& dJn, Real<Prec>& dYn) noexcept {
        using T = Real<Prec>;
        constexpr double xmin = 2;
        constexpr double half = 0.5;
        constexpr double pi_trivial = M_PI;
        const T pi = T(pi_trivial);
        const T epsilon = std::numeric_limits<T>::epsilon();
        const T fpmin = T(std::numeric_limits<T>::min()) / epsilon;
        assert(!n.isNegative() && x.isPositive());

        const Integer nl = [x, n]() noexcept {
            if (x < T(xmin))
                return Integer(T(n + T(half)));
            return Integer(T(n - x + T(xmin - half))).isPositive() ? Integer(T(n - x + T(xmin - half))) : Integer(0);
        }();
        const T mu = n - T(nl);
        assert(mu < x);
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
                sign = sign ^ d.isNegative(); // If d is negative, sign = -sign
                if (abs(T(delta - T(1))) <= epsilon)
                    break;
            } while (true);
        }
        // get f at \nu = \mu
        const T a_Jv_1 = sign ? fpmin : -fpmin; // a is a constant, 1 is the first iteration
        const T a_dJv_1 = f_n * a_Jv_1;
        T a_Jv_m = a_Jv_1;
        T a_dJv_m = a_dJv_1;
        /* Loop from \nu = n to \nu = \mu */ {
            T factor = n * reciprocal_x;
            for (Integer i = nl - 1; !i.isNegative(); --i) {
                T temp = fma(factor, a_Jv_m, a_dJv_m);
                factor -= reciprocal_x;
                a_dJv_m = fma(factor, temp, -a_Jv_m);
                a_Jv_m = std::move(temp);
            }
            if (a_Jv_m.isZero())
                a_Jv_m = epsilon;
        }
        const T f = a_dJv_m / a_Jv_m;

        T Ymu, Ymu_1, Jmu;
        if (x < T(xmin)) { // Temme series method for continued fraction 2
            const T x_2 = x * T(0.5);
            const T pimu = pi * mu;
            const T factor = abs(pimu) < epsilon ? T(1) : pimu / sin(pimu);
            T d = -ln(x_2);
            T e = mu * d;
            const T factor2 = abs(e) < epsilon ? T(1) : sinh(e) / e;
            T gamma1, gamma2, gamma_plus, gamma_minus;
            Internal::besselChebyshevHelper(mu, gamma1, gamma2, gamma_plus, gamma_minus);
            T ff = T(2 / pi_trivial) * factor * fma(gamma1, cosh(e), gamma2 * factor2 * d);
            e = exp(e);
            T p = e / (gamma_plus * pi);
            T q = reciprocal(e * pi * gamma_minus);
            const T pimu_2 = pimu * T(0.5);
            const T factor3 = abs(pimu_2) < epsilon ? T(1) : sin(pimu_2) / pimu_2;
            T r = pi * pimu_2 * square(factor3);
            T c(1);
            d = -square(x_2);
            T sum = fma(r, q, ff);
            T sum1 = p;
            Integer i = 1;
            while (true) {
                const T scalar_i = T(i);
                ff = T(scalar_i * ff + p + q) / (square(scalar_i) - square_mu);
                c *= d / i;
                p /= T(scalar_i - mu);
                q /= T(scalar_i + mu);
                T delta = c * fma(r, q, ff);
                sum += delta;
                sum1 += fma(c, p, -scalar_i * delta);
                if (abs(delta) < (T(1) + abs(sum)) * epsilon)
                    break;
                ++i;
            }
            Ymu = -sum;
            Ymu_1 = -sum1 * reciprocal_x_2;
            const T dYmu = fma(mu * reciprocal_x, Ymu, -Ymu_1);
            Jmu = wronskian / fma(-f, Ymu, dYmu);
        }
        else { // Lentz method for continued fraction 2
            T a = T(0.25) - square_mu;
            T p = reciprocal_x * T(-0.5);
            T q(1);
            const T br = T(2) * x;
            T bi(2);
            // Start from i = 1
            T factor = a * reciprocal_x / fma(p, p, q * q);
            T cr = fma(q, factor, br);
            T ci = fma(p, factor, bi);
            T den = fma(br, br, bi * bi);
            T dr = br / den;
            T di = -bi / den;
            T dlr = fma(cr, dr, -ci * di);
            T dli = fma(cr, di, ci * dr);
            T temp = fma(p, dlr, -q * dli);
            q = fma(p, dli, q * dlr);
            p = temp;
            Integer i = 1;
            while (true) {
                a += T(i << 1U);
                bi += T(2);
                dr = fma(a, dr, br);
                di = fma(a, di, bi);
                if ((abs(dr) + abs(di)) < fpmin)
                    dr = fpmin;
                factor = a / fma(cr, cr, ci * ci);
                cr = fma(cr, factor, br);
                ci = bi - ci * factor;
                if ((abs(cr) + abs(ci)) < fpmin)
                    cr = fpmin;
                den = fma(dr, dr, di * di);
                dr /= den;
                di /= -den;
                dlr = fma(cr, dr, -ci * di);
                dli = fma(cr, di, ci * dr);
                temp = fma(p, dlr, -q * dli);
                q = fma(p, dli, q * dlr);
                p = temp;
                if ((abs(T(dlr - T(1))) + abs(dli)) < epsilon)
                    break;
                ++i;
            }
            const T gamma = (p - f) / q;
            Jmu = sqrt(wronskian / fma(p - f, gamma, q));
            Jmu = a_Jv_m.isPositive() ? Jmu : -Jmu;
            Ymu = Jmu * gamma;
            const T dYmu = Ymu * (p + q / gamma);
            Ymu_1 = fma(mu * reciprocal_x, Ymu, -dYmu);
        }
        const T factor = Jmu / a_Jv_m;
        Jn = a_Jv_1 * factor;
        dJn = a_dJv_1 * factor;
        // Loop from v = \mu to v = n
        T Yv = Ymu, Yv_1 = Ymu_1;
        for (Integer i = 1; i <= nl; ++i) {
            T temp = fma((mu + T(i)) * reciprocal_x_2, Yv_1, -Yv);
            Yv = Yv_1;
            Yv_1 = std::move(temp);
        }
        Yn = Yv;
        dYn = fma(n * reciprocal_x, Yv, -Yv_1);
    }

    template<FloatPrec Prec>
    Real<Prec> besselJn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return Jn;
    }

    template<FloatPrec Prec>
    Real<Prec> besseldJn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return dJn;
    }

    template<FloatPrec Prec>
    Real<Prec> besselYn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return Yn;
    }

    template<FloatPrec Prec>
    Real<Prec> besseldYn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(n, x, Jn, Yn, dJn, dYn);
        return dYn;
    }

    template<FloatPrec Prec>
    void sphericalBesselJn_Yn_dJn_dYn(const Real<Prec>& n, const Real<Prec>& x, Real<Prec>& jn, Real<Prec>& yn, Real<Prec>& djn, Real<Prec>& dyn) noexcept {
        using T = Real<Prec>;
        assert(!n.isNegative() && x.isPositive());

        const T sqrt_pi_2(1.2533141373155002512);
        const T factor = sqrt_pi_2 / sqrt(x);
        const T order = n + T(0.5);
        T Jn, Yn, dJn, dYn;
        besselJn_Yn_dJn_dYn(order, x, Jn, Yn, dJn, dYn);
        jn = factor * Jn;
        yn = factor * Yn;
        djn = factor * (dJn - Jn / (x * T(2)));
        dyn = factor * (dYn - Yn / (x * T(2)));
    }

    template<FloatPrec Prec>
    Real<Prec> sphericalBesselJn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        assert(!n.isNegative() && x.isPositive());

        const T sqrt_pi_2(1.2533141373155002512);
        const T factor = sqrt_pi_2 / sqrt(x);
        const T order = n + T(0.5);
        T Jn, Yn, dJn, dYn;
        besselJn_Yn_dJn_dYn(order, x, Jn, Yn, dJn, dYn);
        return factor * Jn;
    }

    template<FloatPrec Prec>
    Real<Prec> sphericalBesseldJn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        assert(!n.isNegative() && x.isPositive());

        const T sqrt_pi_2(1.2533141373155002512);
        const T factor = sqrt_pi_2 / sqrt(x);
        const T order = n + T(0.5);
        T Jn, Yn, dJn, dYn;
        besselJn_Yn_dJn_dYn(order, x, Jn, Yn, dJn, dYn);
        return factor * (dJn - Jn / (x * T(2)));
    }

    template<FloatPrec Prec>
    Real<Prec> sphericalBesselYn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        assert(!n.isNegative() && x.isPositive());

        const T sqrt_pi_2(1.2533141373155002512);
        const T factor = sqrt_pi_2 / sqrt(x);
        const T order = n + T(0.5);
        T Jn, Yn, dJn, dYn;
        besselJn_Yn_dJn_dYn(order, x, Jn, Yn, dJn, dYn);
        return factor * Yn;
    }

    template<FloatPrec Prec>
    Real<Prec> sphericalBesseldYn(const Real<Prec>& n, const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        assert(!n.isNegative() && x.isPositive());

        const T sqrt_pi_2(1.2533141373155002512);
        const T factor = sqrt_pi_2 / sqrt(x);
        const T order = n + T(0.5);
        T Jn, Yn, dJn, dYn;
        besselJn_Yn_dJn_dYn(order, x, Jn, Yn, dJn, dYn);
        return factor * (dYn - Yn / (x * T(2)));
    }
}
