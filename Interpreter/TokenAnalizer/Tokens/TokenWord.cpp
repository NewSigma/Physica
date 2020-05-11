/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <cstring>
#include "Interpreter/Header/TokenWord.h"

using namespace Physica::Interpreter;

TokenWord::TokenWord(const char* str) : Token(Identifier), data(new char[strlen(str)]) {
    strcpy(data, str);
}

TokenWord::TokenWord(const char* str, int len) : Token(Identifier), data(new char[len]) {
    strcpy(data, str);
}

TokenWord::~TokenWord() {
    delete[] data;
}