/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <iostream>
#include "Physica/Core/IO/Gnuplot.h"
#include "Physica/Core/Exception/BadFileFormatException.h"

using namespace Physica;

std::ostream& Physica::operator<<(std::ostream& os, const Gnuplot& gnuplot) {
    assert(gnuplot.xDatas.getLength() == gnuplot.yDatas.getLength());
    for (size_t i = 0; i < gnuplot.xDatas.getLength(); ++i) {
        auto& xData = gnuplot.xDatas[i];
        auto& yData = gnuplot.yDatas[i];
        assert(xData.getLength() == yData.getLength());
        for (size_t j = 0; j < xData.getLength(); ++j)
            os << xData[j] << ' ' << yData[j] << '\n';
        os << '\n';
    }
    return os;
}

std::istream& Physica::operator>>(std::istream& is, Gnuplot& gnuplot) {
    Array<float64> xBuffer{};
    Array<float64> yBuffer{};
    while (is.good()) {
        int ch = is.peek();
        while (ch == ' ') {
            is.get();
            ch = is.peek();
        }

        if (ch == EOF) {
            if (xBuffer.getLength() != 0) {
                gnuplot.xDatas.append(std::move(xBuffer));
                gnuplot.yDatas.append(std::move(yBuffer));
            }
            break;
        }

        const bool isEmptyLine = ch == '#' || ch == '\n';
        if (isEmptyLine) {
            if (xBuffer.getLength() != 0) {
                gnuplot.xDatas.append(xBuffer);
                gnuplot.yDatas.append(yBuffer);
                xBuffer.resize(0);
                yBuffer.resize(0);
            }
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        double temp1{}, temp2{};
        is >> temp1 >> temp2;
        xBuffer.append(temp1);
        yBuffer.append(temp2);
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (!is)
        throw BadFileFormatException("Gnuplot");
    return is;
}
