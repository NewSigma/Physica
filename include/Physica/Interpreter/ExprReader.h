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
#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include <iosfwd>
#include <Physica/Core/MultiPrecition/Scalar.h>

using Physica::Core::MultiScalar;

namespace Physica::Interpreter {
    class ExprReader {
    private:
        std::list<std::wstring> anti_poland;
    public:
        explicit ExprReader(const std::wstring& str);
        MultiScalar calc();
    private:
        static bool isSign(wchar_t c);
    };
}

#endif
