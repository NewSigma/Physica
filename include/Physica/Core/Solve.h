#ifndef PHYSICA_SOLVE_H
#define PHYSICA_SOLVE_H

namespace Physica::Core {
    class Scalar;

    class Solve {
    public:
        static Scalar bisectionMethod(Scalar func(const Scalar&), const Scalar& n
                , const Scalar& x1, const Scalar& x2);
        static Scalar bisectionMethod(Scalar func(const Scalar&), const Scalar& n
                , const Scalar& x_left, const Scalar& x_right
                , const Scalar& y_left, const Scalar& y2);
    };
}

#endif
