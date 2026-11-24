/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <immintrin.h>

int main() {
    unsigned long long l{};
    return _rdseed64_step(&l);
}
