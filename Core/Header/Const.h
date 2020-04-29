#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

class RealNum;
class Numerical;

class BasicConst {
    static BasicConst* instance;

    const Numerical* plotPoints;
    const Numerical* expectedRelativeError;
    const Numerical* stepSize;
    const Numerical* R_MAX;
    const Numerical* _0;
    const Numerical* _1;
    const Numerical* Minus_1;
    const Numerical* _2;
    const Numerical* Minus_2;
    const Numerical* _3;
    const Numerical* Minus_3;
    const Numerical* _4;
    const Numerical* Minus_4;
public:
    const int GlobalPrecision;
    const int MaxPower;

    BasicConst();
    ~BasicConst();
    static void init();
    static void deInit();
    inline static const BasicConst& getInstance() { return *instance; }

    inline const Numerical& getPlotPoints() const { return *plotPoints; }
    inline const Numerical& getExpectedRelativeError() const { return *expectedRelativeError; }
    inline const Numerical& getStepSize() const { return *stepSize; }
    inline const Numerical& getR_MAX() const { return *R_MAX; }
    inline const Numerical& get_0() const { return *_0; }
    inline const Numerical& get_1() const { return *_1; }
    inline const Numerical& getMinus_1() const { return *Minus_1; }
    inline const Numerical& get_2() const { return *_2; }
    inline const Numerical& getMinus_2() const { return *Minus_2; }
    inline const Numerical& get_3() const { return *_3; }
    inline const Numerical& getMinus_3() const { return *Minus_3; }
    inline const Numerical& get_4() const { return *_4; }
    inline const Numerical& getMinus_4() const { return *Minus_4; }
};

class MathConst {
    static MathConst* instance;

    const RealNum* stepSize;
    const RealNum* _0;
    const RealNum* _1;
    const RealNum* _2;
    const Numerical* PI;
    const Numerical* E;
    //Here PI_2 stands by PI / 2.
    const Numerical* PI_2;
    const Numerical* Minus_PI_2;
public:
    MathConst();
    ~MathConst();
    static void init();
    static void deInit();
    inline static const MathConst& getInstance() { return *instance; }

    inline const RealNum& getStepSize() const { return *stepSize; }
    inline const RealNum& get_0() const { return *_0; }
    inline const RealNum& get_1() const { return *_1; }
    inline const RealNum& get_2() const { return *_2; }
    inline const Numerical& getPI() const { return *PI; }
    inline const Numerical& getE() const { return *E; }
    inline const Numerical& getPI_2() const { return *PI_2; }
    inline const Numerical& getMinus_PI_2() const { return *Minus_PI_2; }
private:
    static Numerical calcPI(int precision);
};

#endif