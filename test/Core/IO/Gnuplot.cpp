/*
 * Copyright 2021-2024 Weibo He.
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
#include <fstream>
#include "Physica/Core/IO/Gnuplot.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

using namespace Physica;

const static char* data1 = "# Test Gnuplot\n"
                           "1.0 1.1\n"
                           "1.1 1.2\n"
                           "\n"
                           "1.0 1.1\n"
                           "1.1 1.2\n"
                           "1.2 1.3\n"
                           "\n"
                           "\n"
                           "0 0\n";

Gnuplot readTest() {
    auto tmp = TempFile("/tmp/tmpXXXXXX");
    std::ofstream os(tmp.getName());
    os << data1;
    os.close();

    Gnuplot gnu{};
    std::ifstream is(tmp.getName());
    is >> gnu;
    is.close();
    return gnu;
}

int main() {
    Gnuplot gnu = readTest();

    const auto& xDatas = gnu.getXDatas();
    if (xDatas.getLength() != 3)
        return 1;
    if (xDatas[0].getLength() != 2 || xDatas[1].getLength() != 3 || xDatas[2].getLength() != 1)
        return 1;
    return 0;
}
