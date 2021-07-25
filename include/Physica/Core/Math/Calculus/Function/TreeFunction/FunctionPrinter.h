/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#ifndef PHYSICA_FUNCTIONPRINTER_H
#define PHYSICA_FUNCTIONPRINTER_H

#include <iosfwd>
#include "TreeFunction.h"

namespace Physica::Core {
    template<ScalarOption option = MultiPrecision, bool errorTrack = true>
    class FunctionPrinter {
        const TreeFunctionData<option, errorTrack>& f;
        std::ostream& os;
    public:
        FunctionPrinter(const TreeFunctionData<option, errorTrack>& f_, std::ostream& os);
        void print() { printImpl(f.getTree()); }
    private:
        void printImpl(const TreeFunctionData<option, errorTrack>& functionTree);
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/FunctionPrinterImpl.h"

#endif