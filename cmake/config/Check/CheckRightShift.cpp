/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
/**
 * Check if this machine uses logical left shift. Logical right shift is unsupported.
 */
int main() {
    return ~0 != (~0 >> 1U);
}
