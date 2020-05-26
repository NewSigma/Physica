/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "FFT.h"

#include <utility>

namespace Physica::Core {
    ///////////////////////////////////FFTBase////////////////////////////////////
    /*
     * Requirements:
     * The elements in y should be sorted.
     *
     * Warning:
     * Length(x) <= Length(y) is required. (Though Length(x) < Length(y) is the often case.)
     * 1 < Length(x) <= Max(SignedNumericalUnit) is required.
     */
    FFTBase::FFTBase(Vector x, Vector y, NormalizationMethod method)
    : data(RowMatrix(new Vector[2]{std::move(x), std::move(y)}, 2)) {
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
    Numerical FFTBase::operator()(const Numerical& n) {
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
    FFT::FFT(Vector x, Vector y, NormalizationMethod method) : FFTBase(std::move(x), std::move(y), method) {
        const auto length = data.column();
        auto& vector_x = data[0];
        /* Handle y */ {
            const auto width = vector_x[length - 1] - vector_x[0];
            const Numerical distance = width / Numerical(static_cast<SignedNumericalUnit>(length - 1));
            data[1] *= distance;
        }
        /* Handle x */
        const Numerical averagePhrase =
                (MathConst::getInstance().getPI() << 1) / Numerical(static_cast<SignedNumericalUnit>(length));
        //There is no need to calculate averagePhrase * 1.
        vector_x[0] = averagePhrase;
        for(size_t i = 1; i < length; ++i)
            vector_x[i] = averagePhrase * Numerical(static_cast<SignedNumericalUnit>(i + 1));
    }
    ///////////////////////////////////InvFFT////////////////////////////////////
    InvFFT::InvFFT(Vector x, Vector y, NormalizationMethod method) : FFTBase(std::move(x), std::move(y), method) {
        const auto length = data.column();
        const Numerical numericalLength(static_cast<SignedNumericalUnit>(length));
        auto& vector_x = data[0];
        /* Handle y */ {
            const auto width = vector_x[length - 1] - vector_x[0];
            const Numerical distance = width / Numerical(static_cast<SignedNumericalUnit>(length - 1));
            data[1] *= distance;
            data[1] /= numericalLength;
        }
        /* Handle x */
        const Numerical averagePhrase =
                -(MathConst::getInstance().getPI() << 1) / numericalLength;
        //There is no need to calculate averagePhrase * 1.
        vector_x[0] = averagePhrase;
        for(size_t i = 1; i < length; ++i)
            vector_x[i] = averagePhrase * Numerical(static_cast<SignedNumericalUnit>(i + 1));
    }
}