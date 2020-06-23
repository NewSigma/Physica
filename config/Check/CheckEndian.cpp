/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */

//Return 0 if little endian, or 1 if big endian.
int main() {
    int i = 1;
    char* p = reinterpret_cast<char*>(&i);
    return *p == 0;
}