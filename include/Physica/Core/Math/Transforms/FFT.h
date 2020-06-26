/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FFT_H
#define PHYSICA_FFT_H

#include "Physica/Core/Math/LinearAlgebra/Matrix/RowMatrix.h"

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
        Scalar operator()(const Scalar& n);
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
