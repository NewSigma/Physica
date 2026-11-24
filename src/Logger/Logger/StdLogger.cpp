/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Logger/Logger/StdLogger.h"
#include <ostream>

using namespace Physica;

StdLogger::StdLogger(std::ostream& stream) : os(stream) {}

void StdLogger::log(LogBuffer& buffer) {
    os << buffer.makeMsgString() << '\n';
}
