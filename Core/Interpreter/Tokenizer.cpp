/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <cstring>
#include <iostream>
#include "../Header/Tokenizer.h"

Tokenizer::Tokenizer(const char* str) : state(Start) {
    int i = 0;
    while(str[i] != '\0') {
        readChar(str[i]);
        ++i;
    }
    resetBuffer(getCurrentType());
}

void Tokenizer::readChar(const char& ch) {
    if(ch == ' ') {
        resetBuffer(getCurrentType());
        return;
    }

    switch(state) {
        case Start:
            switch(ch) {
                case '0':
                    buffer.push_back(ch);
                    state = PreInt;
                    break;
                case '1' ... '9':
                    buffer.push_back(ch);
                    state = Int;
                    break;
                default:
                    std::cout << "TokenizerError" << std::endl;
                    return;
                    //TODO TokenizerError
            }
            break;
        case PreInt:
            switch(ch) {
                case '0':
                    break;
                case '1' ... '9':
                    buffer[0] = ch;
                    state = Int;
                    break;
                case '.':
                    buffer.push_back(ch);
                    state = Float;
                    break;
                default:;
                    resetBuffer(Token::Number);
            }
            break;
        case Int:
            switch(ch) {
                case '0' ... '9':
                    buffer.push_back(ch);
                    break;
                case '.':
                    buffer.push_back(ch);
                    state = Float;
                    break;
                default:;
                    resetBuffer(Token::Number);
            }
            break;
        case Float:
            if(isdigit(ch))
                buffer.push_back(ch);
            else
                resetBuffer(Token::Number);
            break;
    }
}

void Tokenizer::resetBuffer(Token::TokenType type) {
    if(!buffer.empty()) {
        auto temp = new char[buffer.size() + 1];
        std::strcpy(temp, buffer.c_str());
        tokens.emplace_back(temp, type);
        buffer.clear();
        state = Start;
    }
}

Token::TokenType Tokenizer::getCurrentType() {
    switch(state) {
        case Start:
            return Token::Null;
        case PreInt:
        case Int:
        case Float:
            return Token::Number;
    }
}