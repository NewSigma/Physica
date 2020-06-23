/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
int main() {
    int a = -1;
    //Verify two's complement.
    int b = ~a + 1;
    if(b != -a)
        return 1;
    //Verify cast between two's complement and unsigned dose not change the bits.
    unsigned int c = static_cast<unsigned int>(a);
    if((a & c) != a)
        return 1;
    return 0;
}