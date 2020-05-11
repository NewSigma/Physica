#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

namespace Physica::Core {
    class Numerical;

    Numerical rectangular(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1);
    Numerical ladder(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1);
    Numerical simpson(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1);
}

#endif
