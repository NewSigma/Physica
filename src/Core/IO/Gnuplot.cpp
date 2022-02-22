/*
 * Copyright 2021 WeiBo He.
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
#include <iostream>
#include "Physica/Core/IO/Gnuplot.h"
#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica::Core {
    std::ostream& operator<<(std::ostream& os, const Gnuplot& gnuplot) {
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

    std::istream& operator>>(std::istream& is, Gnuplot& gnuplot) {
        typename Gnuplot::VectorType xBuffer{};
        typename Gnuplot::VectorType yBuffer{};
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

            double temp1, temp2;
            is >> temp1 >> temp2;
            xBuffer.append(temp1);
            yBuffer.append(temp2);
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        if (!is)
            throw BadFileFormatException();
        return is;
    }
}
