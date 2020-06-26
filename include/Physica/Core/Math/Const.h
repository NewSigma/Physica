/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CONST_H
#define PHYSICA_CONST_H

#include <iosfwd>
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    constexpr int GlobalPrecision = 4;

    class BasicConst {
        static BasicConst* instance;

        const Scalar<MultiPrecision, false>* plotPoints;
        const Scalar<MultiPrecision, false>* expectedRelativeError;
        const Scalar<MultiPrecision, false>* stepSize;
        const Scalar<MultiPrecision, false>* R_MAX;
        const Scalar<MultiPrecision, false>* _0;
        const Scalar<MultiPrecision, false>* _1;
        const Scalar<MultiPrecision, false>* Minus_1;
        const Scalar<MultiPrecision, false>* _2;
        const Scalar<MultiPrecision, false>* Minus_2;
        const Scalar<MultiPrecision, false>* _3;
        const Scalar<MultiPrecision, false>* Minus_3;
        const Scalar<MultiPrecision, false>* _4;
        const Scalar<MultiPrecision, false>* Minus_4;
        const Scalar<MultiPrecision, false>* _10;
    public:
        const int MaxPower;

        BasicConst();
        ~BasicConst();
        static void init();
        static void deInit();
        inline static const BasicConst& getInstance() { return *instance; }

        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getPlotPoints() const { return *plotPoints; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getExpectedRelativeError() const { return *expectedRelativeError; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getStepSize() const { return *stepSize; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getR_MAX() const { return *R_MAX; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_0() const { return *_0; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_1() const { return *_1; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getMinus_1() const { return *Minus_1; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_2() const { return *_2; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getMinus_2() const { return *Minus_2; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_3() const { return *_3; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getMinus_3() const { return *Minus_3; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_4() const { return *_4; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& getMinus_4() const { return *Minus_4; }
        [[nodiscard]] inline const Scalar<MultiPrecision, false>& get_10() const { return *_10; }
    };

    class MathConst {
        static MathConst* instance;

        const MultiScalar* PI;
        const MultiScalar* E;
        //Here PI_2 stands by PI / 2.
        const MultiScalar* PI_2;
        const MultiScalar* Minus_PI_2;
    public:
        MathConst();
        ~MathConst();
        static void init();
        static void deInit();
        inline static const MathConst& getInstance() { return *instance; }

        [[nodiscard]] inline const MultiScalar& getPI() const { return *PI; }
        [[nodiscard]] inline const MultiScalar& getE() const { return *E; }
        [[nodiscard]] inline const MultiScalar& getPI_2() const { return *PI_2; }
        [[nodiscard]] inline const MultiScalar& getMinus_PI_2() const { return *Minus_PI_2; }
    private:
        static MultiScalar calcPI(int precision);
    };
}

#endif