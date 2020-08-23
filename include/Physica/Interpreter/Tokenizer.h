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
#ifndef PHYSICA_TOKENIZER_H
#define PHYSICA_TOKENIZER_H

#include <list>
#include <string>
#include "Token.h"

namespace Physica::Interpreter {
    class Tokenizer {
        enum NumBufferState {
            NumStart,
            Zero,
            Int,
            Float,
            PreExp,
            SignedPreExp,
            Exp
        };

        enum WordBufferState {
            WordStart,
            Identifier
        };
        std::list<Token*> tokens;
        const char* str;
        //Defined to make debug easier.
        int line;
    public:
        explicit Tokenizer(const char* str);
        ~Tokenizer();
    private:
        void readToken();
        bool readChar(char ch);
        void readNum();
        void readWord();
    };
}

#endif