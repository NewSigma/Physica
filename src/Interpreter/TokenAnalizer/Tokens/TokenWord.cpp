/*
 * Copyright 2020 WeiBo He.
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
#include <cstring>
#include "Physica/Interpreter/TokenWord.h"

using namespace Physica::Interpreter;

TokenWord::TokenWord(const char* str) : Token(Identifier), data(new char[strlen(str)]) {
    strcpy(data, str);
}

TokenWord::TokenWord(const char* str, int len) : Token(Identifier), data(new char[len]) {
    strcpy(data, str);
}

TokenWord::~TokenWord() {
    delete[] data;
}