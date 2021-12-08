/*
 * Copyright 2020-2021 WeiBo He.
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

namespace Physica::Core {
    /*!
     * Improve:
     * 1.Integrate methods may be various.
     * 2.Sampling distance does not have to be equal.
     * 3.A lot of other features.
     */
    template<class ScalarType>
    DFT<ScalarType>::DFT(const Vector<ScalarType>& data_, const ScalarType& distance_)
            : data(data_), distance(distance_) {
        assert(data.getLength() % 2 == 0);
    }

    template<class ScalarType>
    DFT<ScalarType>::DFT(const DFT& dft) : data(dft.data), distance(dft.distance) {}

    template<class ScalarType>
    DFT<ScalarType>::DFT(DFT&& dft) noexcept : data(std::move(dft.data)), distance(std::move(dft.distance)) {}

    template<class ScalarType>
    DFT<ScalarType>& DFT<ScalarType>::operator=(DFT dft) {
        swap(dft);
        return *this;
    }

    template<class ScalarType>
    inline void DFT<ScalarType>::transform() {
        transformImpl(RealType(-2 * M_PI / data.getLength()));
        data *= distance;
    }

    template<class ScalarType>
    inline void DFT<ScalarType>::invTransform() {
        transformImpl(RealType(2 * M_PI / data.getLength()));
        data /= distance * data.getLength();
    }

    template<class ScalarType>
    void DFT<ScalarType>::swap(DFT& dft) {
        swap(data, dft.data);
        swap(distance, dft.distance);
    }

    template<class ScalarType>
    typename DFT<ScalarType>::ComplexType DFT<ScalarType>::getComponent(ssize_t index) const {
        const size_t length = data.getLength();
        assert(length <= SSIZE_MAX);
        assert(index <= static_cast<ssize_t>(length) / 2);
        assert(-static_cast<ssize_t>(length) / 2 <= index);
        if (index < 0)
            index += static_cast<ssize_t>(length);
        return data[index];
    }

    template<class ScalarType>
    Vector<typename DFT<ScalarType>::ComplexType> DFT<ScalarType>::getComponents() const {
        Vector<ComplexType> result = Vector<ComplexType>(data.getLength());
        const ssize_t half_size = static_cast<ssize_t>(data.getLength() / 2);
        for (ssize_t i = -half_size; i < half_size; ++i)
            result[i + half_size] = getComponent(i);
        return result;
    }

    template<class ScalarType>
    void DFT<ScalarType>::transformImpl(const RealType& phase) {
        const size_t length = data.getLength();
        Vector<ComplexType> buffer(length);
        //Optimize:
        //1.i and j is changeable.(dynamic programming)
        //2.Use the formula such as sin(a + b) to avoid calculate sin and cos directly.
        for(size_t i = 0; i < length; ++i) {
            const RealType phase1 = phase * i;
            auto result_i = ComplexType::Zero();
            for(size_t j = 0; j < length; ++j) {
                const RealType phase2 = phase1 * j;
                result_i += ComplexType(cos(phase2), sin(phase2)) * data[j];
            }
            buffer[i] = result_i;
        }
        data = std::move(buffer);
    }
}
