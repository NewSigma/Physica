#include "../../Header/Integrate.h"

Number* Integrate(Number* func(Number*), Number* x0, Number* x1, IntegrateMethod method) {
    switch(method) {
        case Ladder:
            return ladder(func, x0, x1);
        default:
            return rectangular(func, x0, x1);
    }
}

Number* rectangular(Number* func(Number*), Number* x0, Number* x1) {

}

Number* ladder(Number* func(Number*), Number* x0, Number* x1) {

}