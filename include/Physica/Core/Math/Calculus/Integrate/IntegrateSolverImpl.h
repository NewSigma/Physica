/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_INTEGRATESOLVERIMPL_H
#define PHYSICA_INTEGRATESOLVERIMPL_H
/*!
 * Bug: if the start of integrate domain is much larger than step size, the result will be 0. May be use taylor series
 * and expend the function to the first order.
 */
namespace Physica::Core {
    //////////////////////////////////Rectangular//////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    IntegrateSolver<Rectangular, 1, option, errorTrack>::IntegrateSolver(Scalar<option, false> stepSize)
            : stepSize(std::move(stepSize)) {}

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> IntegrateSolver<Rectangular, 1, option, errorTrack>::solve(
            const Integrate<1, option, errorTrack> &i) const {
        const Scalar<option, errorTrack>& from = i.getFrom();
        const Scalar<option, errorTrack>& to = i.getTo();
        const TreeFunction<option, errorTrack>& f = i.getFunction();

        Scalar<option, errorTrack> result = 0;

        Scalar<option, errorTrack> start(from);
        while(start < to) {
            result += f(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Ladder//////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    IntegrateSolver<Ladder, 1, option, errorTrack>::IntegrateSolver(Scalar<option, false> stepSize)
            : stepSize(std::move(stepSize)) {}

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> IntegrateSolver<Ladder, 1, option, errorTrack>::solve(
            const Integrate<1, option, errorTrack> &i) const {
        const Scalar<option, errorTrack>& from = i.getFrom();
        const Scalar<option, errorTrack>& to = i.getTo();
        const TreeFunctionData<option, errorTrack>& f = i.getFunction();

        Scalar<option, errorTrack> result = ((f(from) + f(to)) >> 1);
        Scalar<option, errorTrack> start(from + stepSize);
        while(start < to) {
            result += f(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Simpson//////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    IntegrateSolver<Simpson, 1, option, errorTrack>::IntegrateSolver(Scalar<option, false> stepSize)
            : stepSize(std::move(stepSize)) {}

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> IntegrateSolver<Simpson, 1, option, errorTrack>::solve(
            const Integrate<1, option, errorTrack> &i) const {
        const Scalar<option, errorTrack>& from = i.getFrom();
        const Scalar<option, errorTrack>& to = i.getTo();
        const TreeFunctionData<option, errorTrack>& f = i.getFunction();

        Scalar<option, errorTrack> result = f(from) + f(to);
        const auto& _0 = BasicConst::getInstance()._0;
        Scalar<option, errorTrack> odd(_0);
        Scalar<option, errorTrack> even(_0);
        bool b = true;
        Scalar<option, errorTrack> start = from + stepSize;
        while(start < to) {
            Scalar<option, errorTrack>& toChange = b ? odd : even;
            b = !b;
            toChange += f(start);
            start += stepSize;
        }
        odd <<= 2;
        even <<= 1;
        result += odd + even;
        result *= stepSize;
        result /= BasicConst::getInstance()._3;
        return result;
    }
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    IntegrateSolver<Tanh_Sinh, 1, option, errorTrack>::IntegrateSolver(Scalar<option, false> stepSize, size_t pointCount)
            : stepSize(std::move(stepSize)), pointCount(pointCount) {}
    /*!
     * Reference:
     * [1] Vanherck, Joren Sorée, Bart Magnus, Wim.
     * Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic. http://arxiv.org/abs/2007.15057
     */
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> IntegrateSolver<Tanh_Sinh, 1, option, errorTrack>::solve(
            const Integrate<1, option, errorTrack>& i) const {
        const Scalar<option, errorTrack>& from = i.getFrom();
        const Scalar<option, errorTrack>& to = i.getTo();

        const Scalar<option, errorTrack> constant1 = (to - from) >> 1;
        const Scalar<option, errorTrack> constant2 = constant1 + from;
        const auto& PI_2 = MathConst::getInstance().PI_2;

        const TreeFunctionData<option, errorTrack>& f = i.getFunction();
        //Integrate value on the origin.
        Scalar<option, errorTrack> result = PI_2 * f(constant2);
        Scalar<option, errorTrack> point_x = 0;
        for(size_t j = 0; j < pointCount; ++j) {
            point_x += stepSize;
            const Scalar<option, errorTrack> PI_2_sinh = PI_2 * sinh(point_x);
            const Scalar<option, errorTrack> cosh_PI_2_sinh = cosh(PI_2_sinh);
            const Scalar<option, errorTrack> phi = sinh(PI_2_sinh) / cosh_PI_2_sinh;
            const Scalar<option, errorTrack> phi_derivative = PI_2 * cosh(point_x) / square(cosh_PI_2_sinh);
            result += phi_derivative * (f(constant2 + constant1 * phi) + f(constant2 - constant1 * phi));
        }
        result *= constant1 * stepSize;
        return result;
    }
}

#endif
