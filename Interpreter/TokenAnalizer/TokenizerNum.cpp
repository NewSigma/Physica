/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include "Interpreter/Header/TokenizerNum.h"
/*
 * Number part of Tokenizer.
 */
//Standard output: (+/-) 12.75 E (+/-) 13
//FIXME Can not recognize +1 or -2 and so on.
TokenizerNum::TokenizerNum(const char*& str) : state(Start) {
    while(readChar(*str))
        ++str;
    if(!stoppable)
        qFatal("Encountered tokenizer error.");
}

bool TokenizerNum::readChar(char ch) {
    switch(state) {
        case Start:
            type = Token::Integer;
            switch(ch) {
                case '0':
                    buffer.push_back(ch);
                    state = Zero;
                    break;
                case '1' ... '9':
                    buffer.push_back(ch);
                    state = Int;
                    break;
                default:
                    return stoppable = false;
            }
            break;
        case Zero:
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
                    type = Token::Float;
                    break;
                case 'E':
                case 'e':
                    buffer.push_back('E');
                    state = PreExp;
                    type = Token::IntegerExp;
                    break;
                default:
                    return false;
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
                    type = Token::Float;
                    break;
                case 'E':
                case 'e':
                    buffer.push_back('E');
                    state = PreExp;
                    type = Token::IntegerExp;
                    break;
                default:
                    return false;
            }
            break;
        case Float:
            switch(ch) {
                case '0' ... '9':
                    buffer.push_back(ch);
                    break;
                case 'E':
                case 'e':
                    buffer.push_back('E');
                    state = PreExp;
                    type = Token::FloatExp;
                    break;
                default:
                    return false;
            }
            break;
        case PreExp:
            stoppable = false;
            switch(ch) {
                case '0':
                    buffer.push_back('+');
                    buffer.push_back('0');
                    state = SignedPreExp;
                    break;
                case '1' ... '9':
                    buffer.push_back('+');
                    buffer.push_back(ch);
                    state = Exp;
                    break;
                case '+':
                case '-':
                    buffer.push_back(ch);
                    state = SignedPreExp;
                    break;
                default:
                    return false;
            }
            break;
        case SignedPreExp:
            stoppable = false;
            switch(ch) {
                case '0':
                    break;
                case '1' ... '9':
                    buffer.push_back(ch);
                    state = Exp;
                    break;
                default:
                    return false;
            }
            break;
        case Exp:
            stoppable = true;
            if(isdigit(ch))
                buffer.push_back(ch);
            else
                return false;
            break;
    }
    return true;
}