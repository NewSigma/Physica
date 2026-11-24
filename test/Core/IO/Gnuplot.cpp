/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <fstream>
#include "Physica/Core/IO/Gnuplot.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;

namespace {
    const char* data1 = "# Test Gnuplot\n"
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
}

int main() {
    Gnuplot gnu = readTest();

    const auto& xDatas = gnu.getXDatas();
    expect(xDatas.getLength() == 3);
    expect(xDatas[0].getLength() == 2);
    expect(xDatas[1].getLength() == 3);
    expect(xDatas[2].getLength() == 1);
    return 0;
}
