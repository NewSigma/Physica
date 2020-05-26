/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FFT_H
#define PHYSICA_FFT_H

#include "Matrix.h"

namespace Physica::Core {
    ///////////////////////////////////FFTBase////////////////////////////////////
    class FFTBase {
    protected:
        RowMatrix data;
    public:
        //NormalizationMethod decides the discretization method used for the integration.
        enum NormalizationMethod {
            RectangularMethod,
            LadderMethod
        };
        FFTBase(FFTBase&& base) noexcept;
        Numerical operator()(const Numerical& n);
    protected:
        FFTBase(Vector x, Vector y, NormalizationMethod method);
    private:
        void ladderMethod();
    };
    ///////////////////////////////////FFT////////////////////////////////////
    class FFT : public FFTBase {
    public:
        FFT(Vector x, Vector y, NormalizationMethod method = RectangularMethod);
    };
    ///////////////////////////////////InvFFT////////////////////////////////////
    class InvFFT : public FFTBase {
    public:
        InvFFT(Vector x, Vector y, NormalizationMethod method = RectangularMethod);
    };
}

#endif
