/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKEN_H
#define PHYSICA_TOKEN_H

#include <string>

class Token {
    char* data;
public:
    enum TokenType {
        Null,
        Integer,
        IntegerExp,
        Float,
        FloatExp,
        Identifier,
        //Operators
        OperatorAddSub,     // +, -
        OperatorMulDiv,     // *, /
        OperatorOpeAssign,  // +=, -=, *=, /= and so on.
        OperatorLess,       // <
        OperatorLessEq,     // <
        OperatorLarger,     // >
        OperatorLargerEq,   // >
        OperatorAssign,     // =
        OperatorEq,
        OperatorNot,
        OperatorNotEq,
        //Keywords
        KeyWordIf,
        KeyWordElse,
        KeyWordSwitch,
        KeyWordFor,
        KeyWordDo,
        KeyWordWhile
        //Boundarys
    } const type;

    Token(const std::string& data, TokenType type);
    ~Token();
};

#endif
