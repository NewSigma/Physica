/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <filesystem>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Exception/IOException.h"

using namespace Physica;

namespace {
    constexpr const char* BadDataFileMsg = "[Error]: Bad data file";
    constexpr const char* BadLabelFileMsg = "[Error]: Bad label file";
}

Mnist::Mnist(const char* folder) {
    std::filesystem::path p(folder);
    p.append("train-images.idx3-ubyte");
    trainSamples = readDatas(p.string());
    if (trainSamples.getLength() != 60000)
        throw BadFileFormatException(BadDataFileMsg);
    p.remove_filename();

    p.append("t10k-images.idx3-ubyte");
    testSamples = readDatas(p.string());
    if (testSamples.getLength() != 10000)
        throw BadFileFormatException(BadDataFileMsg);
    p.remove_filename();

    p.append("train-labels.idx1-ubyte");
    trainLabels = readLabels(p.string());
    if (trainLabels.getLength() != 60000)
        throw BadFileFormatException(BadLabelFileMsg);
    p.remove_filename();

    p.append("t10k-labels.idx1-ubyte");
    testLabels = readLabels(p.string());
    if (testLabels.getLength() != 10000)
        throw BadFileFormatException(BadLabelFileMsg);
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
        throw BadFileFormatException(BadLabelFileMsg);

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
        throw BadFileFormatException(BadLabelFileMsg);

    LabelArray result(readInt(fin));
    fin.read(reinterpret_cast<char*>(result.data()), result.getLength());
    return result;
}

int32_t Mnist::readInt(std::ifstream& fin) {
    IntDecomp temp{};
    fin.read(temp.c.data(), sizeof(temp.c));
    std::swap(temp.c[0], temp.c[3]);
    std::swap(temp.c[1], temp.c[2]);
    return temp.i;
}
