/*
 * Copyright 2024 WeiBo He.
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

#include "SimpleDataset.h"

namespace Physica::Core {
    template<class SampleType, class LabelType>
    class device_obj<SimpleDataset<SampleType, LabelType>> {
        using host_obj = SimpleDataset<SampleType, LabelType>;
        using This = device_obj<host_obj>;
    public:
        using SampleArray = Utils::device_obj<Utils::Array<SampleType>>;
        using LabelArray = Utils::device_obj<Utils::Array<LabelType>>;
        using DataType = std::pair<typename SampleArray::ValueType, typename LabelArray::ValueType>;
    private:
        SampleArray samples;
        LabelArray labels;
    public:
        device_obj() = default;
        device_obj(const host_obj& obj);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] host_obj toHost() const;
        void toHost(host_obj& obj) const;
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __device__ SampleArray& getSamples() noexcept { return samples; }
        [[nodiscard]] __device__ const SampleArray& getSamples() const noexcept { return samples; }
        [[nodiscard]] __device__ LabelArray& getLabels() noexcept { return labels; }
        [[nodiscard]] __device__ const LabelArray& getLabels() const noexcept { return labels; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return samples.getLength(); }
    };

    template<class SampleType, class LabelType>
    device_obj<SimpleDataset<SampleType, LabelType>>::device_obj(const host_obj& obj) : samples(obj.samples), labels(obj.labels) {}

    template<class SampleType, class LabelType>
    SimpleDataset<SampleType, LabelType> device_obj<SimpleDataset<SampleType, LabelType>>::toHost() const {
        SimpleDataset<SampleType, LabelType> result{};
        result.samples = samples.toHost();
        result.labels = labels.toHost();
        return result;
    }

    template<class SampleType, class LabelType>
    void device_obj<SimpleDataset<SampleType, LabelType>>::toHost(host_obj& obj) const {
        samples.toHost(obj.samples);
        labels.toHost(obj.labels);
    }

    template<class SampleType, class LabelType>
    void device_obj<SimpleDataset<SampleType, LabelType>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        samples.swap(obj.samples);
        labels.swap(obj.labels);
    }

    template<class SampleType, class LabelType>
    typename SimpleDataset<SampleType, LabelType>::device_obj_type SimpleDataset<SampleType, LabelType>::toDevice() const {
        return device_obj_type(*this);
    }

    template<class SampleType, class LabelType>
    void SimpleDataset<SampleType, LabelType>::toDevice(device_obj_type& obj) const {
        samples.toDevice(obj.samples);
        labels.toDevice(obj.labels);
    }
}
