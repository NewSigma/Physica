/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <immintrin.h>

int main() {
    unsigned long long l{};
    return _rdrand64_step(&l);
}
