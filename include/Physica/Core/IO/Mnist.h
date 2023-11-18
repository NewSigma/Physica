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

        DataArray trainSamples;
        DataArray testSamples;
        LabelArray trainLabels;
        LabelArray testLabels;
    public:
        Mnist(const char* folder);
        Mnist(const Mnist&) = default;
        Mnist(Mnist&&) noexcept = default;
        ~Mnist() = default;
        /* Operators */
        Mnist& operator=(Mnist obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class ScalarType> [[nodiscard]] inline Vector<ScalarType> makeTrainSample(size_t i) const;
        template<class ScalarType> [[nodiscard]] inline Vector<ScalarType> makeTestSample(size_t i) const;
        template<class ScalarType> [[nodiscard]] inline Vector<ScalarType> makeTrainLabel(size_t i) const;
        template<class ScalarType> [[nodiscard]] inline Vector<ScalarType> makeTestLabel(size_t i) const;
        void swap(Mnist& obj) noexcept;
        /* Getters */
        [[nodiscard]] const DataArray& getTrainSamples() const noexcept { return trainSamples; }
        [[nodiscard]] const DataArray& getTestSamples() const noexcept { return testSamples; }
        [[nodiscard]] const LabelArray& getTrainLabels() const noexcept { return trainLabels; }
        [[nodiscard]] const LabelArray& getTestLabels() const noexcept { return testLabels; }
        [[nodiscard]] size_t getNumTrainSample() const noexcept { return trainSamples.getLength(); }
        [[nodiscard]] size_t getNumTestSample() const noexcept { return testSamples.getLength(); }
        /* Static members */
        template<class ScalarType>
        [[nodiscard]] static Vector<ScalarType, NumCategory> makeLabelVector(unsigned char label);
    private:
        [[nodiscard]] static DataArray readDatas(const std::string& path);
        [[nodiscard]] static LabelArray readLabels(const std::string& path);
        [[nodiscard]] static int32_t readInt(std::ifstream& fin);
    };

    template<class ScalarType>
    inline Vector<ScalarType> Mnist::makeTrainSample(size_t i) const {
        assert(i < getNumTrainSample() && "[Error]: Index overflow");
        return getTrainSamples()[i].asVector<ScalarType>();
    }

    template<class ScalarType>
    inline Vector<ScalarType> Mnist::makeTestSample(size_t i) const {
        assert(i < getNumTestSample() && "[Error]: Index overflow");
        return getTestSamples()[i].asVector<ScalarType>();
    }

    template<class ScalarType>
    inline Vector<ScalarType> Mnist::makeTrainLabel(size_t i) const {
        return makeLabelVector<ScalarType>(getTrainLabels()[i]);
    }

    template<class ScalarType>
    inline Vector<ScalarType> Mnist::makeTestLabel(size_t i) const {
        return makeLabelVector<ScalarType>(getTestLabels()[i]);
    }

    template<class ScalarType>
    Vector<ScalarType> Mnist::ImageType::asVector() const {
        Vector<ScalarType> result(NumPixelInImage);
        for (size_t i = 0; i < NumPixelInImage; ++i)
            result[i] = ScalarType(pixels[i]);
        return result;
    }

    template<class ScalarType>
    Vector<ScalarType, Mnist::NumCategory> Mnist::makeLabelVector(unsigned char label) {
        using PlainScalar = typename ScalarType::PlainScalar;
        assert(label < NumCategory && "[Error]: Invalid category");
        Vector<PlainScalar, Mnist::NumCategory> result(NumCategory, 0);
        result[label] = PlainScalar(1);
        return result;
    }
}
