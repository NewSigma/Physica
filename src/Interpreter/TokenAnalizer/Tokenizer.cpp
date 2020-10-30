/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <QtCore/qlogging.h>
#include <Physica/Core/MultiPrecition/Scalar.h>
#include "Physica/Interpreter/TokenWord.h"
#include "Physica/Interpreter/TokenNum.h"
#include "Physica/Interpreter/Tokenizer.h"
#include "Physica/Core/MultiPrecition/Const.h"

using namespace Physica::Interpreter;
using namespace Physica::Core;

Tokenizer::Tokenizer(const char* str) : str(str), line(1) {
    Token::init();
    while(*str != '\0') {
        Q_UNUSED(str)
        readToken();
    }
}

Tokenizer::~Tokenizer() {
    for(auto p : tokens) {
        switch(p->type) {
            case Token::Numeric:
            case Token::Identifier:
                delete p;
            default:;
        }
    }
    Token::deInit();
}

void Tokenizer::readToken() {
    const char ch = *str;

    switch(ch) {
        case '\n':
            ++line;
        case ' ':
        case '\t':
            ++str;
            return;
        //Handle operators and boundary symbols
        case '+':
            if(readChar('='))
                tokens.push_back(Token::operatorAddEq);
            else
                tokens.push_back(Token::operatorAdd);
            return;
        case '-':
            if(readChar('='))
                tokens.push_back(Token::operatorSubEq);
            else
                tokens.push_back(Token::operatorSub);
            return;
        case '/': {
            //Handle notes
            const char next = *(str + 1);
            if(next == '/') {
                str += 2;
                while(*str != '\n')
                    ++str;
                return;
            }
            if(next == '*') {
                str += 2;
                while(*str != '*' || *(str + 1) != '/')
                    ++str;
                str += 2;
                return;
            }
        }
        case '*':
            if(readChar('='))
                tokens.push_back(Token::operatorMulEq);
            else
                tokens.push_back(Token::operatorMul);
            return;
        case '<':
            if(readChar('='))
                tokens.push_back(Token::operatorLessEq);
            else
                tokens.push_back(Token::operatorLess);
            return;
        case '>':
            if(readChar('='))
                tokens.push_back(Token::operatorLargerEq);
            else
                tokens.push_back(Token::operatorLarger);
            return;
        case '=':
            if(readChar('='))
                tokens.push_back(Token::operatorEq);
            else
                tokens.push_back(Token::operatorAssign);
            return;
        case '!':
            if(readChar('='))
                tokens.push_back(Token::operatorNotEq);
            else
                tokens.push_back(Token::operatorNot);
            return;
        default:;
    }

    if(isdigit(ch)) {
        readNum();
        return;
    }

    if(isalpha(ch) || ch == '_') {
        readWord();
        return;
    }
    qFatal("Unrecognized token.");
}
//Return if the next char equals to ch.
bool Tokenizer::readChar(char ch) {
    bool result = *(str + 1) == ch;
    if(result)
        ++str;
    return result;
}

void Tokenizer::readNum() {
    bool stoppable = true, go_on = true;
    NumBufferState state = NumStart;
    Scalar<MultiPrecision, false> num(BasicConst::getInstance()._0);
    Scalar exp(BasicConst::getInstance()._0);
    SignedScalarUnit power = 0;
    bool isExpPositive = true;
    while(go_on) {
        switch(state) {
            case NumStart:
                switch(*str) {
                    case '0':
                        state = Zero;
                        ++str;
                        break;
                    case '1' ... '9':
                        num = Scalar<MultiPrecision, false>(*str - '0');
                        state = Int;
                        ++str;
                        break;
                    default:
                        go_on = stoppable = false;
                }
                break;
            case Zero:
                switch(*str) {
                    case '0':
                        break;
                    case '1' ... '9':
                        num = Scalar<MultiPrecision, false>(*str - '0');
                        state = Int;
                        ++str;
                        break;
                    case '.':
                        state = Float;
                        break;
                    case 'E':
                    case 'e':
                        state = PreExp;
                        break;
                    default:
                        go_on = false;
                }
                break;
            case Int:
                switch(*str) {
                    case '0' ... '9':
                        num *= BasicConst::getInstance()._10;
                        num += MultiScalar(static_cast<SignedScalarUnit>(*str - '0'));
                        ++str;
                        break;
                    case '.':
                        state = Float;
                        break;
                    case 'E':
                    case 'e':
                        state = PreExp;
                        break;
                    default:
                        go_on = false;
                }
                break;
            case Float:
                switch(*str) {
                    case '0' ... '9':
                        num *= BasicConst::getInstance()._10;
                        num += MultiScalar(static_cast<SignedScalarUnit>(*str - '0'));
                        --power;
                        ++str;
                        break;
                    case 'E':
                    case 'e':
                        state = PreExp;
                        break;
                    default:
                        go_on = false;
                }
                break;
            case PreExp:
                stoppable = false;
                switch(*str) {
                    case '+':
                    case '0':
                        state = SignedPreExp;
                        ++str;
                        break;
                    case '1' ... '9':
                        exp += MultiScalar(static_cast<SignedScalarUnit>(*str - '0'));
                        state = Exp;
                        ++str;
                        break;
                    case '-':
                        isExpPositive = false;
                        state = SignedPreExp;
                        ++str;
                        break;
                    default:
                        go_on = false;
                }
                break;
            case SignedPreExp:
                stoppable = false;
                switch(*str) {
                    case '0':
                        break;
                    case '1' ... '9':
                        exp *= BasicConst::getInstance()._10;
                        exp += MultiScalar(static_cast<SignedScalarUnit>(*str - '0'));
                        state = Exp;
                        ++str;
                        break;
                    default:
                        go_on = false;
                }
                break;
            case Exp:
                stoppable = true;
                if(isdigit(*str)) {
                    exp *= BasicConst::getInstance()._10;
                    exp += MultiScalar(static_cast<SignedScalarUnit>(*str - '0'));
                    ++str;
                }
                else
                    go_on = false;
                break;
        }
    }
    if(!stoppable)
        qFatal("Encountered tokenizer error.");
    if(!isExpPositive)
        exp.toOpposite();
    num *= (BasicConst::getInstance()._10 ^ (exp + MultiScalar(power)));
    tokens.push_back(new TokenNum(num));
}

void Tokenizer::readWord() {
    bool stoppable = true, go_on = true;
    WordBufferState state = WordBufferState::WordStart;
    std::string buffer{};
    while(go_on) {
        switch(state) {
            case WordStart:
                if(isalpha(*str) || *str == '_') {
                    buffer.push_back(*str);
                    state = Identifier;
                    ++str;
                }
                else
                    go_on = stoppable = false;
                break;
            case Identifier:
                if(isalpha(*str) || isdigit(*str) || *str == '_') {
                    buffer.push_back(*str);
                    state = Identifier;
                    ++str;
                }
                else
                    go_on = false;
                break;
        }
    }
    if(!stoppable)
        qFatal("Encountered tokenizer error.");
    //Recognize keywords
    if(buffer == "if")
        tokens.push_back(Token::keyWordIf);
    else if(buffer == "else")
        tokens.push_back(Token::keyWordElse);
    else if(buffer == "switch")
        tokens.push_back(Token::keyWordSwitch);
    else if(buffer == "for")
        tokens.push_back(Token::keyWordFor);
    else if(buffer == "do")
        tokens.push_back(Token::keyWordDo);
    else if(buffer == "while")
        tokens.push_back(Token::keyWordWhile);
    else
        tokens.push_back(new TokenWord(buffer.c_str(), static_cast<int>(buffer.length()) + 1));
}