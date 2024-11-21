/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/Utils/DirStack.h"

namespace Physica::Core {
    Mnist::Mnist(const char* folder) {
        DirStack dir(folder);
        dir.push("train-images.idx3-ubyte");
        trainSamples = readDatas(dir.toPath());
        if (trainSamples.getLength() != 60000)
            throw BadFileFormatException("[Error]: Bad data file");
        dir.pop();

        dir.push("t10k-images.idx3-ubyte");
        testSamples = readDatas(dir.toPath());
        if (testSamples.getLength() != 10000)
            throw BadFileFormatException("[Error]: Bad data file");
        dir.pop();

        dir.push("train-labels.idx1-ubyte");
        trainLabels = readLabels(dir.toPath());
        if (trainLabels.getLength() != 60000)
            throw BadFileFormatException("[Error]: Bad label file");
        dir.pop();

        dir.push("t10k-labels.idx1-ubyte");
        testLabels = readLabels(dir.toPath());
        if (testLabels.getLength() != 10000)
            throw BadFileFormatException("[Error]: Bad label file");
    }

    void Mnist::swap(Mnist& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        trainSamples.swap(obj.trainSamples);
        testSamples.swap(obj.testSamples);
        trainLabels.swap(obj.trainLabels);
        testLabels.swap(obj.testLabels);
    }

    Mnist::DataArray Mnist::readDatas(const std::string& path) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: File not found");
        const auto magic = readInt(fin);
        const auto numData = readInt(fin);
        const auto row = readInt(fin);
        const auto col = readInt(fin);
        if (magic != 2051 || row != ImageSize || col != ImageSize)
            throw BadFileFormatException("[Error]: Bad label file");
        
        DataArray result(numData);
        fin.read(reinterpret_cast<char*>(result.data()), result.getLength() * sizeof(ImageType));
        return result;
    }

    Mnist::LabelArray Mnist::readLabels(const std::string& path) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: File not found");
        const auto magic = readInt(fin);
        if (magic != 2049)
            throw BadFileFormatException("[Error]: Bad label file");
        
        LabelArray result(readInt(fin));
        fin.read(reinterpret_cast<char*>(result.data()), result.getLength());
        return result;
    }

    int32_t Mnist::readInt(std::ifstream& fin) {
        IntDecomp temp{};
        fin.read(temp.c, sizeof(temp.c));
        std::swap(temp.c[0], temp.c[3]);
        std::swap(temp.c[1], temp.c[2]);
        return temp.i;
    }
}
