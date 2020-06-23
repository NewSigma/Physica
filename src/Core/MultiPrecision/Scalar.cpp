/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    void Scalar<0, true>::swap(Scalar<0, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<0, false>::swap(s);
    }

    void Scalar<1, true>::swap(Scalar<1, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<1, false>::swap(s);
    }
}