#include "../../Header/Integrate.h"

AbstractNum* Integrate(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1, IntegrateMethod method) {
    switch(method) {
        case Ladder:
            return ladder(func, x0, x1);
        default:
            return rectangular(func, x0, x1);
    }
}

AbstractNum* rectangular(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1) {
    return nullptr;
}

AbstractNum* ladder(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1) {
    return nullptr;
}