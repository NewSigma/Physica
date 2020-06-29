/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BISECTIONMETHOD_H
#define PHYSICA_BISECTIONMETHOD_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2) {
        Scalar<MultiPrecision, false> y1 = func(x1);
        Scalar<MultiPrecision, false> y2 = func(x2);
        return bisectionMethod<errorTrack>(func, n, x1, x2, y1, y2);
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2
            , const Scalar<MultiPrecision, false>& y1, const Scalar<MultiPrecision, false>& y2);

    template<>
    Scalar<MultiPrecision, false> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2
            , const Scalar<MultiPrecision, false>& y1, const Scalar<MultiPrecision, false>& y2);

    template<>
    Scalar<MultiPrecision, true> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2
            , const Scalar<MultiPrecision, false>& y1, const Scalar<MultiPrecision, false>& y2);
}

#endif