#ifndef PHYSICA_SOLVE_H
#define PHYSICA_SOLVE_H

namespace Physica::Core {
    class Numerical;

    class Solve {
    public:
        static Numerical bisectionMethod(Numerical func(const Numerical&), const Numerical& n
                , const Numerical& x1, const Numerical& x2);
        static Numerical bisectionMethod(Numerical func(const Numerical&), const Numerical& n
                , const Numerical& x_left, const Numerical& x_right
                , const Numerical& y_left, const Numerical& y2);
    };
}

#endif
