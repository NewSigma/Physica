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
#ifndef PHYSICA_TREEFUNCTIONPRINTER_H
#define PHYSICA_TREEFUNCTIONPRINTER_H

#include <list>
#include <iosfwd>
#include "TreeFunctionData.h"

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    class TreeFunction;

    namespace {
        class TreeFunctionPrinterPrivate {
        protected:
            static const char* rightBranch;
            static const char* leftBranch;
            static const char* straightLine;
            static const char* space;
            static const char* spaces;
        };
    }

    template<ScalarType type, bool errorTrack>
    class TreeFunctionPrinter : private TreeFunctionPrinterPrivate {
        std::list<const char*> list;
        const TreeFunction<type, errorTrack>& f;
        std::ostream& os;
    public:
        TreeFunctionPrinter(const TreeFunction<type, errorTrack>& f, std::ostream& os);
        void print() { printImpl(f.getTree(), false); } //true or false leads to the same result.
    private:
        void printImpl(const TreeFunctionData<type, errorTrack>& functionTree, bool isLeft);
        void printList();
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionDataPrinterImpl.h"

#endif