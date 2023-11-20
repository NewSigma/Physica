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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    template<class SampleType, class LabelType>
    class SimpleDataset {
        using ScalarType1 = typename SampleType::ScalarType;
        using ScalarType2 = typename LabelType::ScalarType;
        static_assert(std::is_same<ScalarType1, ScalarType2>::value, "[Error]: ScalarType should match");
        static_assert(!ScalarType1::isDifferentiable, "[Error]: Data in a dataset must not be differentiable");
        using SampleArray = Utils::Array<SampleType>;
        using LabelArray = Utils::Array<LabelType>;
    private:
        Utils::Array<SampleType> samples;
        Utils::Array<LabelType> labels;
    public:
        SimpleDataset(SampleArray samples_, LabelArray labels_);
        SimpleDataset(const SimpleDataset&) = default;
        SimpleDataset(SimpleDataset&&) noexcept = default;
        ~SimpleDataset() = default;
        /* Operators */
        SimpleDataset& operator=(SimpleDataset obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(SimpleDataset& obj) noexcept;
        /* Getters */
        [[nodiscard]] const SampleArray& getSamples() const noexcept { return samples; }
        [[nodiscard]] const LabelArray& getLabels() const noexcept { return labels; }
        [[nodiscard]] size_t getSize() const noexcept { return samples.getLength(); }
    };

    template<class SampleType, class LabelType>
    SimpleDataset<SampleType, LabelType>::SimpleDataset(SampleArray samples_, LabelArray labels_)
            : samples(std::move(samples_)), labels(std::move(labels_)) {
        assert(samples.getLength() == labels.getLength() && "[Error]: Length of samples and labels do not match");
    }

    template<class SampleType, class LabelType>
    void SimpleDataset<SampleType, LabelType>::swap(SimpleDataset& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        labels.swap(obj.labels);
    }
}
