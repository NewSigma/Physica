/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include <Interpreter/Header/TokenizerNum.h>
#include <Interpreter/Header/TokenizerWord.h>
#include "Interpreter/Header/Tokenizer.h"

Tokenizer::Tokenizer(const char* str) : str(str), line(1) {
    while(*str != '\0') {
        Q_UNUSED(str)
        readWord();
    }
}

void Tokenizer::readWord() {
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
        case '-':
            handleOperator(Token::OperatorAddSub, Token::OperatorOpeAssign);
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
            handleOperator(Token::OperatorMulDiv, Token::OperatorOpeAssign);
            return;
        case '<':
            handleOperator(Token::OperatorLess, Token::OperatorLessEq);
            return;
        case '>':
            handleOperator(Token::OperatorLarger, Token::OperatorLargerEq);
            return;
        case '=':
            handleOperator(Token::OperatorAssign, Token::OperatorEq);
            return;
        case '!':
            handleOperator(Token::OperatorNot, Token::OperatorNotEq);
            return;
        default:;
    }

    if(isdigit(ch)) {
        TokenizerNum tokenizerNum(str);
        tokens.emplace_back(tokenizerNum.getBuffer(), tokenizerNum.getType());
        return;
    }

    if(isalpha(ch) || ch == '_') {
        TokenizerWord tokenizerId(str);
        tokens.emplace_back(tokenizerId.getBuffer(), tokenizerId.getType());
        return;
    }

    qFatal("Unrecognized token.");
}
//If the operator is accompanied with '=', use pair. Or use single.
void Tokenizer::handleOperator(Token::TokenType single, Token::TokenType pair) {
    std::string buffer{*str};
    ++str;
    if(*str == '=') {
        buffer.push_back('=');
        ++str;
        tokens.emplace_back(buffer, pair);
        return;
    }
    tokens.emplace_back(buffer, single);
}