/*
 * Copyright 2023 WeiBo He.
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
#pragma once

#include <string>
#include <fstream>
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] http://yann.lecun.com/exdb/mnist/
     */
    class Mnist {
        constexpr static size_t ImageSize = 28;
        struct ImageType {
            unsigned char pixels[ImageSize * ImageSize];
        };

        union IntDecomp {
            char c[4];
            int32_t i;
        };

        using LabelArray = Utils::Array<unsigned char>;
        using DataArray = Utils::Array<ImageType>;
    private:
        LabelArray trainLabels;
        LabelArray testLabels;
        DataArray trainDatas;
        DataArray testDatas;
    public:
        Mnist(const char* folder);
        Mnist(const Mnist&) = default;
        Mnist(Mnist&&) noexcept = default;
        ~Mnist() = default;
        /* Operators */
        Mnist& operator=(Mnist obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(Mnist& obj) noexcept;
        /* Getters */
        [[nodiscard]] const LabelArray& getTrainLabels() const noexcept { return trainLabels; }
        [[nodiscard]] const LabelArray& getTestLabels() const noexcept { return testLabels; }
        [[nodiscard]] size_t getNumTrainLabel() const noexcept { return trainLabels.getLength(); }
        [[nodiscard]] size_t getNumTestLabel() const noexcept { return testLabels.getLength(); }
        [[nodiscard]] const DataArray& getTrainDatas() const noexcept { return trainDatas; }
        [[nodiscard]] const DataArray& getTestDatas() const noexcept { return testDatas; }
    private:
        [[nodiscard]] static LabelArray readLabels(const std::string& path);
        [[nodiscard]] static DataArray readDatas(const std::string& path);
        [[nodiscard]] static int32_t readInt(std::ifstream& fin);
    };
}
