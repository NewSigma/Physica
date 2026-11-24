/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<Inverse<M>>> : public Traits<Inverse<M>> {};
}
