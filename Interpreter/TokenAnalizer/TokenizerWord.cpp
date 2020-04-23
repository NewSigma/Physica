/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
/*
 * Identifier part of Tokenizer.
 */
#include <QtCore/qlogging.h>
#include "TokenizerWord.h"

TokenizerWord::TokenizerWord(const char*& str) : state(Start) {
    while(readChar(*str))
        ++str;
    if(!stoppable)
        qFatal("Encountered tokenizer error.");
    //Recognize keywords
    if(buffer == "if")
        type = Token::KeyWordIf;
    else if(buffer == "else")
        type = Token::KeyWordElse;
    else if(buffer == "switch")
        type = Token::KeyWordSwitch;
    else if(buffer == "for")
        type = Token::KeyWordFor;
    else if(buffer == "do")
        type = Token::KeyWordDo;
    else if(buffer == "while")
        type = Token::KeyWordWhile;
}

bool TokenizerWord::readChar(char ch) {
    switch(state) {
        case Start:
            type = Token::Identifier;
            if(isalpha(ch) || ch == '_') {
                buffer.push_back(ch);
                state = Identifier;
            }
            else
                return stoppable = false;
            break;
        case Identifier:
            if(isalpha(ch) || isdigit(ch) || ch == '_') {
                buffer.push_back(ch);
                state = Identifier;
            }
            else
                return false;
            break;
    }
    return true;
}