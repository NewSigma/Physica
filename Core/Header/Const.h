#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

namespace Physica::Core {
    class RealNum;
    class Scalar;

    class BasicConst {
        static BasicConst* instance;

        const Scalar* plotPoints;
        const Scalar* expectedRelativeError;
        const Scalar* stepSize;
        const Scalar* R_MAX;
        const Scalar* _0;
        const Scalar* _1;
        const Scalar* Minus_1;
        const Scalar* _2;
        const Scalar* Minus_2;
        const Scalar* _3;
        const Scalar* Minus_3;
        const Scalar* _4;
        const Scalar* Minus_4;
        const Scalar* _10;
    public:
        const int GlobalPrecision;
        const int MaxPower;

        BasicConst();
        ~BasicConst();
        static void init();
        static void deInit();
        inline static const BasicConst& getInstance() { return *instance; }

        [[nodiscard]] inline const Scalar& getPlotPoints() const { return *plotPoints; }
        [[nodiscard]] inline const Scalar& getExpectedRelativeError() const { return *expectedRelativeError; }
        [[nodiscard]] inline const Scalar& getStepSize() const { return *stepSize; }
        [[nodiscard]] inline const Scalar& getR_MAX() const { return *R_MAX; }
        [[nodiscard]] inline const Scalar& get_0() const { return *_0; }
        [[nodiscard]] inline const Scalar& get_1() const { return *_1; }
        [[nodiscard]] inline const Scalar& getMinus_1() const { return *Minus_1; }
        [[nodiscard]] inline const Scalar& get_2() const { return *_2; }
        [[nodiscard]] inline const Scalar& getMinus_2() const { return *Minus_2; }
        [[nodiscard]] inline const Scalar& get_3() const { return *_3; }
        [[nodiscard]] inline const Scalar& getMinus_3() const { return *Minus_3; }
        [[nodiscard]] inline const Scalar& get_4() const { return *_4; }
        [[nodiscard]] inline const Scalar& getMinus_4() const { return *Minus_4; }
        [[nodiscard]] inline const Scalar& get_10() const { return *_10; }
    };

    class MathConst {
        static MathConst* instance;

        const RealNum* stepSize;
        const RealNum* _0;
        const RealNum* _1;
        const RealNum* _2;
        const Scalar* PI;
        const Scalar* E;
        //Here PI_2 stands by PI / 2.
        const Scalar* PI_2;
        const Scalar* Minus_PI_2;
    public:
        MathConst();
        ~MathConst();
        static void init();
        static void deInit();
        inline static const MathConst& getInstance() { return *instance; }

        [[nodiscard]] inline const RealNum& getStepSize() const { return *stepSize; }
        [[nodiscard]] inline const RealNum& get_0() const { return *_0; }
        [[nodiscard]] inline const RealNum& get_1() const { return *_1; }
        [[nodiscard]] inline const RealNum& get_2() const { return *_2; }
        [[nodiscard]] inline const Scalar& getPI() const { return *PI; }
        [[nodiscard]] inline const Scalar& getE() const { return *E; }
        [[nodiscard]] inline const Scalar& getPI_2() const { return *PI_2; }
        [[nodiscard]] inline const Scalar& getMinus_PI_2() const { return *Minus_PI_2; }
    private:
        static Scalar calcPI(int precision);
    };
}

#endif