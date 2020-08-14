/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_INTEGRATESOLVERIMPL_H
#define PHYSICA_INTEGRATESOLVERIMPL_H
/*!
 * Bug: if the start of integrate domain is much larger than step size, the result will be 0. May be use taylor series
 * and expend the function to the first order.
 */
namespace Physica::Core {
    //////////////////////////////////Rectangular//////////////////////////////////
    template<ScalarType type, bool errorTrack>
    IntegrateSolver<Rectangular, 1, type, errorTrack>::IntegrateSolver(Scalar<type, false> stepSize) : stepSize(std::move(stepSize)) {}

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> IntegrateSolver<Rectangular, 1, type, errorTrack>::solve(
            const Integrate<1, type, errorTrack> &i) const {
        const Scalar<type, errorTrack>& from = i.getFrom();
        const Scalar<type, errorTrack>& to = i.getTo();
        const TreeFunction<type, errorTrack>& f = i.getFunction();

        Scalar<type, errorTrack> result(BasicConst::getInstance().get_0());

        Scalar<type, errorTrack> start(from);
        while(start < to) {
            result += func(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Ladder//////////////////////////////////
    template<ScalarType type, bool errorTrack>
    IntegrateSolver<Ladder, 1, type, errorTrack>::IntegrateSolver(Scalar<type, false> stepSize) : stepSize(std::move(stepSize)) {}

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> IntegrateSolver<Ladder, 1, type, errorTrack>::solve(
            const Integrate<1, type, errorTrack> &i) const {
        const Scalar<type, errorTrack>& from = i.getFrom();
        const Scalar<type, errorTrack>& to = i.getTo();
        const TreeFunction<type, errorTrack>& f = i.getFunction();

        Scalar<type, errorTrack> result = ((f(from) + f(to)) >> 1);
        Scalar<type, errorTrack> start(from + stepSize);
        while(start < to) {
            result += f(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Simpson//////////////////////////////////
    template<ScalarType type, bool errorTrack>
    IntegrateSolver<Simpson, 1, type, errorTrack>::IntegrateSolver(Scalar<type, false> stepSize) : stepSize(std::move(stepSize)) {}

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> IntegrateSolver<Simpson, 1, type, errorTrack>::solve(
            const Integrate<1, type, errorTrack> &i) const {
        const Scalar<type, errorTrack>& from = i.getFrom();
        const Scalar<type, errorTrack>& to = i.getTo();
        const TreeFunction<type, errorTrack>& f = i.getFunction();

        Scalar<type, errorTrack> result = f(from) + f(to);
        Scalar<type, errorTrack> odd(BasicConst::getInstance().get_0());
        Scalar<type, errorTrack> even(BasicConst::getInstance().get_0());
        bool b = true;
        Scalar<type, errorTrack> start = from + stepSize;
        while(start < to) {
            Scalar<type, errorTrack>& toChange = b ? odd : even;
            b = !b;
            toChange += f(start);
            start += stepSize;
        }
        odd <<= 2;
        even <<= 1;
        result += odd + even;
        result *= stepSize;
        result /= BasicConst::getInstance().get_3();
        return result;
    }
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    template<ScalarType type, bool errorTrack>
    IntegrateSolver<Tanh_Sinh, 1, type, errorTrack>::IntegrateSolver(Scalar<type, false> stepSize, size_t pointCount)
            : stepSize(std::move(stepSize)), pointCount(pointCount) {}
    /*!
     * Reference:
     * [1] Vanherck, Joren Sorée, Bart Magnus, Wim.
     * Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic. http://arxiv.org/abs/2007.15057
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> IntegrateSolver<Tanh_Sinh, 1, type, errorTrack>::solve(
            const Integrate<1, type, errorTrack>& i) const {
        const Scalar<type, errorTrack>& from = i.getFrom();
        const Scalar<type, errorTrack>& to = i.getTo();

        const Scalar<type, errorTrack> constant1 = (to - from) >> 1;
        const Scalar<type, errorTrack> constant2 = constant1 + from;
        const auto& PI_2 = MathConst::getInstance().getPI_2();

        const TreeFunction<type, errorTrack>& f = i.getFunction();
        //Integrate value on the origin.
        Scalar<type, errorTrack> result = PI_2 * f(constant2);
        Scalar<type, errorTrack> point_x = 0;
        for(size_t j = 0; j < pointCount; ++j) {
            point_x += stepSize;
            const Scalar<type, errorTrack> PI_2_sinh = PI_2 * sinh(point_x);
            const Scalar<type, errorTrack> cosh_PI_2_sinh = cosh(PI_2_sinh);
            const Scalar<type, errorTrack> phi = sinh(PI_2_sinh) / cosh_PI_2_sinh;
            const Scalar<type, errorTrack> phi_derivative = PI_2 * cosh(point_x) / square(cosh_PI_2_sinh);
            result += phi_derivative * (f(constant2 + constant1 * phi) + f(constant2 - constant1 * phi));
        }
        result *= constant1 * stepSize;
        return result;
    }
}

#endif
