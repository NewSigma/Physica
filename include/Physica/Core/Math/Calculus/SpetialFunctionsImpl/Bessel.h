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
        const T eight = T(8);
        if (ax < eight) {
            y = square(ax);
            ans1 = T(57568490574.0) + y * (T(-13362590354.0) + y * (T(651619640.7) + y * (T(-11214424.18) + y * (T(77392.33017) + y * T(-184.9052456)))));
            ans2 = T(57568490411.0) + y * (T(1029532985.0) + y * (T(9494680.718) + y * (T(59272.64853) + y * (T(267.8532712) + y))));
            return ans1 / ans2;
        }
        else {
            const T z = eight / ax;
            y = square(z);
            ans1 = T(1) + y * (T(-0.1098628627E-2) + y * (T(0.2734510407E-4) + y * (T(-0.2073370639E-5) + y * T(0.2093887211E-6))));
            ans2 = T(-0.1562499995E-1) + y * (T(0.1430488765E-3) + y * (T(-0.6911147651E-5) + y * (T(0.7621095161E-6) - y * T(0.934945152E-7))));
            const T xx = ax - T(0.785398164);
            return sqrt(T(0.636619772) / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
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
        const T eight = T(8);
        if (ax < eight) {
            y = square(ax);
            ans1 = x * (T(72362614232.0) + y * (T(-7895059235.0) + y * (T(242396853.1) + y * (T(-2972611.439) + y * (T(15704.48260) + y * T(-30.16036606))))));
            ans2 = T(144725228442.0) + y * (T(2300535178.0) + y * (T(18583304.74) + y * (T(99447.43394) + y * (T(376.9991397) + y))));
            return ans1 / ans2;
        }
        else {
            const T z = eight / ax;
            y = square(z);
            ans1 = T(1) + y * (T(0.183105E-2) + y * (T(-0.3516396496E-4) + y * (T(0.2457520174E-5) + y * T(-0.240337019E-6))));
            ans2 = T(0.04687499995) + y * (T(-0.2002690873E-3) + y * (T(0.8449199096E-5) + y * (T(-0.88228987E-6) + y * T(0.105787412E-6))));
            const T xx = ax - T(2.356194491);
            const T result = sqrt(T(0.636619772) / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
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
}