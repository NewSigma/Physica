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
#pragma once

#include "Physica/Core/Utils/Container/Array.h"

namespace Physica::Core {
    template<class SampleType, class LabelType>
    class SimpleDataset {
        using This = SimpleDataset<SampleType, LabelType>;
        using SplitResultType = std::pair<This, This>;
    public:
        using device_obj_type = device_obj<This>;
        using SampleArray = Array<SampleType>;
        using LabelArray = Array<LabelType>;
        using DataType = std::pair<SampleType, LabelType>;
    private:
        SampleArray samples;
        LabelArray labels;
    public:
        SimpleDataset() = default;
        SimpleDataset(SampleArray samples_, LabelArray labels_);
        template<class OtherSample, class OtherLabel>
        SimpleDataset(const SimpleDataset<OtherSample, OtherLabel>& other);
        SimpleDataset(const SimpleDataset&) = default;
        SimpleDataset(SimpleDataset&&) noexcept = default;
        ~SimpleDataset() = default;
        /* Operators */
        SimpleDataset& operator=(SimpleDataset obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] inline DataType operator[](size_t index) const;
        /* Operations */
        inline void reserve(size_t size);
        inline void append(DataType data);
        template<RandomGenerator R>
        SplitResultType randomSplit(size_t firstSize) const;
        void swap(SimpleDataset& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] SampleArray& getSamples() noexcept { return samples; }
        [[nodiscard]] const SampleArray& getSamples() const noexcept { return samples; }
        [[nodiscard]] LabelArray& getLabels() noexcept { return labels; }
        [[nodiscard]] const LabelArray& getLabels() const noexcept { return labels; }
        [[nodiscard]] size_t getSize() const noexcept { return samples.getLength(); }
    private:
        friend class device_obj<This>;
    };

    template<class SampleType, class LabelType>
    SimpleDataset<SampleType, LabelType>::SimpleDataset(SampleArray samples_, LabelArray labels_)
            : samples(std::move(samples_)), labels(std::move(labels_)) {
        assert(samples.getLength() == labels.getLength() && "[Error]: Length of samples and labels do not match");
    }

    template<class SampleType, class LabelType>
    template<class OtherSample, class OtherLabel>
    SimpleDataset<SampleType, LabelType>::SimpleDataset(const SimpleDataset<OtherSample, OtherLabel>& other) {
        const size_t size = other.getSize();
        reserve(size);
        for (size_t i = 0; i < size; ++i) {
            auto pair = other[i];
            append(std::make_pair(OtherSample(std::move(pair.first)), OtherLabel(std::move(pair.second))));
        }
    }

    template<class SampleType, class LabelType>
    inline SimpleDataset<SampleType, LabelType>::DataType
    SimpleDataset<SampleType, LabelType>::operator[](size_t index) const {
        assert(index < getSize() && "[Error]: Index overflow");
        return std::make_pair(samples[index], labels[index]);
    }

    template<class SampleType, class LabelType>
    inline void SimpleDataset<SampleType, LabelType>::reserve(size_t size) {
        samples.reserve(size);
        labels.reserve(size);
    }

    template<class SampleType, class LabelType>
    inline void SimpleDataset<SampleType, LabelType>::append(DataType data) {
        samples.append(std::move(data.first));
        labels.append(std::move(data.second));
    }

    template<class SampleType, class LabelType>
    template<RandomGenerator R>
    SimpleDataset<SampleType, LabelType>::SplitResultType
    SimpleDataset<SampleType, LabelType>::randomSplit(size_t firstSize) const {
        assert(firstSize > 0 && "[Error]: Spliting a zero size dataset does nothing");
        assert(firstSize < getSize() && "[Error]: Split a dataset whose size is larger than original");
        const size_t secondSize = getSize() - firstSize;
        const bool isFirstLarger = firstSize >= secondSize;
        const size_t largeSize = isFirstLarger ? firstSize : secondSize;
        const size_t smallSize = isFirstLarger ? secondSize : firstSize;
        This large{}, small{};
        large.reserve(largeSize);
        small.reserve(smallSize);

        Array<size_t> permutation(getSize());
        for (size_t i = 0; i < getSize(); ++i)
            permutation[i] = i;

        for (size_t i = 0; i < smallSize; ++i) {
            std::uniform_int_distribution<size_t> dist(i, getSize() - 1);
            const size_t j = dist(R::getInstance());
            std::swap(permutation[i], permutation[j]);
            small.append((*this)[j]);
        }
        for (size_t i = smallSize; i < getSize(); ++i)
            large.append((*this)[permutation[i]]);

        if (isFirstLarger)
            return std::make_pair(std::move(large), std::move(small));
        return std::make_pair(std::move(small), std::move(large));
    }

    template<class SampleType, class LabelType>
    void SimpleDataset<SampleType, LabelType>::swap(SimpleDataset& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        labels.swap(obj.labels);
    }
}
