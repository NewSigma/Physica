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
#include "Physica/Interpreter/Token.h"

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