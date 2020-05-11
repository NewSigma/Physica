/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Interpreter/Header/Token.h"

using namespace Physica::Interpreter;

Token* Token::operatorAdd = nullptr; Token* Token::operatorSub = nullptr; Token* Token::operatorMul = nullptr; Token* Token::operatorDiv = nullptr;
Token* Token::operatorAddEq = nullptr; Token* Token::operatorSubEq = nullptr; Token* Token::operatorMulEq = nullptr; Token* Token::operatorDivEq = nullptr;
Token* Token::operatorLess = nullptr; Token* Token::operatorLessEq = nullptr; Token* Token::operatorLarger = nullptr; Token* Token::operatorLargerEq = nullptr;
Token* Token::operatorAssign = nullptr; Token* Token::operatorEq = nullptr; Token* Token::operatorNot = nullptr; Token* Token::operatorNotEq = nullptr;
Token* Token::keyWordIf = nullptr; Token* Token::keyWordElse = nullptr; Token* Token::keyWordSwitch = nullptr;
Token* Token::keyWordFor = nullptr; Token* Token::keyWordDo = nullptr; Token* Token::keyWordWhile = nullptr;

Token::Token(TokenType t) : type(t) {}

void Token::init() {
    operatorAdd = new Token(OperatorAdd);
    operatorSub = new Token(OperatorSub);
    operatorMul = new Token(OperatorMul);
    operatorDiv = new Token(OperatorDiv);
    operatorAddEq = new Token(OperatorAddEq);
    operatorSubEq = new Token(OperatorSubEq);
    operatorMulEq = new Token(OperatorMulEq);
    operatorDivEq = new Token(OperatorDivEq);
    operatorLess = new Token(OperatorLess);
    operatorLessEq = new Token(OperatorLessEq);
    operatorLarger = new Token(OperatorLarger);
    operatorLargerEq = new Token(OperatorLargerEq);
    operatorAssign = new Token(OperatorAssign);
    operatorEq = new Token(OperatorEq);
    operatorNot = new Token(OperatorNot);
    operatorNotEq = new Token(OperatorNotEq);
    keyWordIf = new Token(KeyWordIf);
    keyWordElse = new Token(KeyWordElse);
    keyWordSwitch = new Token(KeyWordSwitch);
    keyWordFor = new Token(KeyWordFor);
    keyWordDo = new Token(KeyWordDo);
    keyWordWhile = new Token(KeyWordWhile);
}

void Token::deInit() {
    delete operatorAdd; delete operatorSub; delete operatorMul; delete operatorDiv;
    delete operatorAddEq; delete operatorSubEq; delete operatorMulEq; delete operatorDivEq;
    delete operatorLess; delete operatorLessEq; delete operatorLarger; delete operatorLargerEq;
    delete operatorAssign; delete operatorEq; delete operatorNot; delete operatorNotEq;
    delete keyWordIf; delete keyWordElse; delete keyWordSwitch;
    delete keyWordFor; delete keyWordDo; delete keyWordWhile;
}