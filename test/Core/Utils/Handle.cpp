/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Utils/Handle.h"

using namespace Physica;

static_assert(std::is_trivially_copyable_v<Handle<HandleType::MPI_Comm>>);

int main() {}
