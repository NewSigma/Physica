/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/ML/NeuralNetwork/SimpleDataset.h"

namespace Physica {
    /**
     * Reference:
     * [1] http://yann.lecun.com/exdb/mnist/
     */
    class PHYSICA_API Mnist {
        using This = Mnist;
    public:
        constexpr static size_t ImageSize = 28;
        constexpr static size_t NumPixelInImage = ImageSize * ImageSize;
        constexpr static size_t NumCategory = 10;
        struct ImageType {
            Array<unsigned char, NumPixelInImage> pixels;
            /* Operations */
            template<Vector V>
            [[nodiscard]] V asVector() const;
        };

        template<class T>
        using DatasetType = SimpleDataset<T, unsigned char>;
    private:
        union IntDecomp {
            Array<char, 4> c;
            int32_t i;
        };

        using LabelArray = Array<unsigned char>;
        using DataArray = Array<ImageType>;

        DataArray trainSamples;
        DataArray testSamples;
        LabelArray trainLabels;
        LabelArray testLabels;
    public:
        Mnist(const char* folder = ".");
        Mnist(const Mnist&) = default;
        Mnist(Mnist&&) noexcept = default;
        ~Mnist() = default;
        /* Operators */
        Mnist& operator=(Mnist obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Vector V> DatasetType<V> makeTrainDataset() const;
        template<Vector V> DatasetType<V> makeTestDataset() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getTrainSamples() const noexcept { return trainSamples; }
        [[nodiscard]] const auto& getTestSamples() const noexcept { return testSamples; }
        [[nodiscard]] const auto& getTrainLabels() const noexcept { return trainLabels; }
        [[nodiscard]] const auto& getTestLabels() const noexcept { return testLabels; }
        [[nodiscard]] size_t getNumTrainSample() const noexcept { return trainSamples.getLength(); }
        [[nodiscard]] size_t getNumTestSample() const noexcept { return testSamples.getLength(); }
    private:
        [[nodiscard]] static DataArray readDatas(const std::string& path);
        [[nodiscard]] static LabelArray readLabels(const std::string& path);
        [[nodiscard]] static int32_t readInt(std::ifstream& fin);
    };

    template<Vector V>
    Mnist::DatasetType<V> Mnist::makeTrainDataset() const {
        using SampleArray = DatasetType<V>::SampleArray;
        const size_t numSample = getNumTrainSample();
        SampleArray samples(numSample);
        for (size_t i = 0; i < numSample; ++i)
            samples[i] = trainSamples[i].asVector<V>();
        return DatasetType<V>(std::move(samples), trainLabels);
    }

    template<Vector V>
    Mnist::DatasetType<V> Mnist::makeTestDataset() const {
        using SampleArray = DatasetType<V>::SampleArray;
        const size_t numSample = getNumTestSample();
        SampleArray samples(numSample);
        for (size_t i = 0; i < numSample; ++i)
            samples[i] = testSamples[i].asVector<V>();
        return DatasetType<V>(std::move(samples), testLabels);
    }

    template<Vector V>
    V Mnist::ImageType::asVector() const {
        using ScalarType = V::ScalarType::ValueType;
        VectorND<ScalarType> result(NumPixelInImage);
        for (size_t i = 0; i < NumPixelInImage; ++i)
            result[i] = ScalarType(pixels[i]);
        return result;
    }
}
