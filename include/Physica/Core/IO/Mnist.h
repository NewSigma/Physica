/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/AI/NeuralNetwork/SimpleDataset.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] http://yann.lecun.com/exdb/mnist/
     */
    class PHYSICA_API Mnist {
    public:
        constexpr static size_t ImageSize = 28;
        constexpr static size_t NumPixelInImage = ImageSize * ImageSize;
        constexpr static size_t NumCategory = 10;
        struct ImageType {
            unsigned char pixels[NumPixelInImage];
            /* Operations */
            template<class VectorType>
            [[nodiscard]] VectorType asVector() const;
        };
        template<class VectorType> using DatasetType = SimpleDataset<VectorType, unsigned char>;
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
        template<class VectorType> DatasetType<VectorType> makeTrainDataset() const;
        template<class VectorType> DatasetType<VectorType> makeTestDataset() const;
        void swap(Mnist& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DataArray& getTrainSamples() const noexcept { return trainSamples; }
        [[nodiscard]] const DataArray& getTestSamples() const noexcept { return testSamples; }
        [[nodiscard]] const LabelArray& getTrainLabels() const noexcept { return trainLabels; }
        [[nodiscard]] const LabelArray& getTestLabels() const noexcept { return testLabels; }
        [[nodiscard]] size_t getNumTrainSample() const noexcept { return trainSamples.getLength(); }
        [[nodiscard]] size_t getNumTestSample() const noexcept { return testSamples.getLength(); }
    private:
        [[nodiscard]] static DataArray readDatas(const std::string& path);
        [[nodiscard]] static LabelArray readLabels(const std::string& path);
        [[nodiscard]] static int32_t readInt(std::ifstream& fin);
    };

    template<class VectorType>
    typename Mnist::DatasetType<VectorType> Mnist::makeTrainDataset() const {
        static_assert(is_vector<VectorType>::value, "[Error]: This is not a vector");
        using SampleArray = typename DatasetType<VectorType>::SampleArray;
        const size_t numSample = getNumTrainSample();
        SampleArray samples(numSample);
        for (size_t i = 0; i < numSample; ++i)
            samples[i] = trainSamples[i].asVector<VectorType>();
        return DatasetType<VectorType>(std::move(samples), trainLabels);
    }

    template<class VectorType>
    typename Mnist::DatasetType<VectorType> Mnist::makeTestDataset() const {
        static_assert(is_vector<VectorType>::value, "[Error]: This is not a vector");
        using SampleArray = typename DatasetType<VectorType>::SampleArray;
        const size_t numSample = getNumTestSample();
        SampleArray samples(numSample);
        for (size_t i = 0; i < numSample; ++i)
            samples[i] = testSamples[i].asVector<VectorType>();
        return DatasetType<VectorType>(std::move(samples), testLabels);
    }

    template<class VectorType>
    VectorType Mnist::ImageType::asVector() const {
        static_assert(is_vector<VectorType>::value, "[Error]: This is not a vector");
        using ScalarType = typename VectorType::ScalarType::PlainScalar;
        Vector<ScalarType> result(NumPixelInImage);
        for (size_t i = 0; i < NumPixelInImage; ++i)
            result[i] = ScalarType(pixels[i]);
        return result;
    }
}
