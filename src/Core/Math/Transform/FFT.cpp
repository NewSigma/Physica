/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Transform/FFTImpl/ThreadGuardFFTW.h"

namespace Physica::Internal {
    ThreadGuardFFTW& ThreadGuardFFTW::getInstance() noexcept {
        constinit static ThreadGuardFFTW instance{};
        return instance;
    }
}
