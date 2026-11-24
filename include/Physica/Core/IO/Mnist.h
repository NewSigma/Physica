/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
