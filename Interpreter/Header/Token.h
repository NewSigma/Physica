/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKEN_H
#define PHYSICA_TOKEN_H

#include <string>

class Token {
public:
    enum TokenType {
        Numeric,
        Identifier,
        //Operators
        OperatorAdd, OperatorSub, OperatorMul, OperatorDiv,
        OperatorAddEq, OperatorSubEq, OperatorMulEq, OperatorDivEq,
        OperatorLess, OperatorLessEq, OperatorLarger, OperatorLargerEq,   // <, <=, >, >=
        OperatorAssign, OperatorEq, OperatorNot, OperatorNotEq, // =, ==, !, !=
        //Keywords
        KeyWordIf, KeyWordElse, KeyWordSwitch,
        KeyWordFor, KeyWordDo, KeyWordWhile
        //Boundarys
    } const type;

    explicit Token(TokenType type);

    static void init();
    static void deInit();

    static Token* operatorAdd, *operatorSub, *operatorMul, *operatorDiv
                , *operatorAddEq, *operatorSubEq, *operatorMulEq, *operatorDivEq
                , *operatorLess, *operatorLessEq, *operatorLarger, *operatorLargerEq
                , *operatorAssign, *operatorEq, *operatorNot, *operatorNotEq
                , *keyWordIf, *keyWordElse, *keyWordSwitch
                , *keyWordFor, *keyWordDo, *keyWordWhile;
};

#endif
