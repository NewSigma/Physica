/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <utility>
#include "Physica/Core/Math/Transforms/FFT.h"

namespace Physica::Core {
    ///////////////////////////////////FFTBase////////////////////////////////////
    /*
     * Requirements:
     * The elements in y should be sorted.
     *
     * Warning:
     * Length(x) <= Length(y) is required. (Though Length(x) < Length(y) is the often case.)
     * 1 < Length(x) <= Max(SignedScalarUnit) is required.
     */
    FFTBase::FFTBase(Vector<MultiScalar> x, Vector<MultiScalar> y, NormalizationMethod method)
    : data(RowMatrix(2, 2)) {
        data.allocate(std::move(x), 0);
        data.allocate(std::move(y), 1);
        //Apply integration.
        switch(method) {
            case LadderMethod:
                ladderMethod();
                break;
            default: //RectangularMethod
                break;
        }
    }

    FFTBase::FFTBase(FFTBase&& base) noexcept : data(std::move(base.data)) {}
    //Calculate value.
    MultiScalar FFTBase::operator()(const MultiScalar& n) {
        const auto real = cos(data[0] * n);
        //Complex number not implemented.
        //const auto imagine = sin(data[0] * n);
        return data[1] * real;
    }

    void FFTBase::ladderMethod() {
        auto& vector_y = data[1];
        const auto lastIndex = vector_y.getLength() - 1;
        for(size_t i = 0; i < lastIndex; ++i) {
            vector_y[i] += vector_y[i + 1];
            vector_y[i] >>= 1;
        }
    }
    ///////////////////////////////////FFT////////////////////////////////////
    FFT::FFT(Vector<MultiScalar> x, Vector<MultiScalar> y, NormalizationMethod method) : FFTBase(std::move(x), std::move(y), method) {
        const auto length = data.column();
        auto& vector_x = data[0];
        /* Handle y */ {
            const auto width = vector_x[length - 1] - vector_x[0];
            const MultiScalar distance = width / MultiScalar(static_cast<SignedScalarUnit>(length - 1));
            data[1] *= distance;
        }
        /* Handle x */
        const MultiScalar averagePhrase =
                (MathConst::getInstance().getPI() << 1) / MultiScalar(static_cast<SignedScalarUnit>(length));
        //There is no need to calculate averagePhrase * 1.
        vector_x[0] = averagePhrase;
        for(size_t i = 1; i < length; ++i)
            vector_x[i] = averagePhrase * MultiScalar(static_cast<SignedScalarUnit>(i + 1));
    }
    ///////////////////////////////////InvFFT////////////////////////////////////
    InvFFT::InvFFT(Vector<MultiScalar> x, Vector<MultiScalar> y, NormalizationMethod method) : FFTBase(std::move(x), std::move(y), method) {
        const auto length = data.column();
        const MultiScalar numericalLength(static_cast<SignedScalarUnit>(length));
        auto& vector_x = data[0];
        /* Handle y */ {
            const auto width = vector_x[length - 1] - vector_x[0];
            const MultiScalar distance = width / MultiScalar(static_cast<SignedScalarUnit>(length - 1));
            data[1] *= distance;
            data[1] /= numericalLength;
        }
        /* Handle x */
        const MultiScalar averagePhrase =
                -(MathConst::getInstance().getPI() << 1) / numericalLength;
        //There is no need to calculate averagePhrase * 1.
        vector_x[0] = averagePhrase;
        for(size_t i = 1; i < length; ++i)
            vector_x[i] = averagePhrase * MultiScalar(static_cast<SignedScalarUnit>(i + 1));
    }
}