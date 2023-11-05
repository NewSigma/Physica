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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] http://yann.lecun.com/exdb/mnist/
     */
    class Mnist {
    public:
        constexpr static size_t ImageSize = 28;
        constexpr static size_t NumPixelInImage = ImageSize * ImageSize;
        constexpr static size_t NumCategory = 10;
        struct ImageType {
            unsigned char pixels[NumPixelInImage];
            /* Operations */
            template<class ScalarType>
            [[nodiscard]] Vector<ScalarType> asVector() const;
        };
    private:
        union IntDecomp {
            char c[4];
            int32_t i;
        };

        using LabelArray = Utils::Array<unsigned char>;
        using DataArray = Utils::Array<ImageType>;

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
        [[nodiscard]] const DataArray& getTrainDatas() const noexcept { return trainDatas; }
        [[nodiscard]] const DataArray& getTestDatas() const noexcept { return testDatas; }
        [[nodiscard]] size_t getNumTrainData() const noexcept { return trainDatas.getLength(); }
        [[nodiscard]] size_t getNumTestData() const noexcept { return testDatas.getLength(); }
        /* Static members */
        template<class ScalarType>
        [[nodiscard]] static Vector<ScalarType, NumCategory> makeLabelVector(unsigned char label);
    private:
        [[nodiscard]] static LabelArray readLabels(const std::string& path);
        [[nodiscard]] static DataArray readDatas(const std::string& path);
        [[nodiscard]] static int32_t readInt(std::ifstream& fin);
    };

    template<class ScalarType>
    Vector<ScalarType> Mnist::ImageType::asVector() const {
        Vector<ScalarType> result(NumPixelInImage);
        for (size_t i = 0; i < NumPixelInImage; ++i)
            result[i] = ScalarType(pixels[i]);
        return result;
    }

    template<class ScalarType>
    Vector<ScalarType, Mnist::NumCategory> Mnist::makeLabelVector(unsigned char label) {
        assert(label < NumCategory && "[Error]: Invalid category");
        Vector<ScalarType, Mnist::NumCategory> result(NumCategory, 0);
        result[label] = ScalarType(1);
        return result;
    }
}
