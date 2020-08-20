/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
/*!
 * Check if this machine uses logical left shift. Logical right shift is unsupported.
 */
int main() {
    return ~0 != (~0 >> 1U);
}