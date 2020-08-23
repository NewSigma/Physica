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
#ifndef PHYSICA_TOKEN_H
#define PHYSICA_TOKEN_H

#include <string>

namespace Physica::Interpreter {
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
}

#endif
