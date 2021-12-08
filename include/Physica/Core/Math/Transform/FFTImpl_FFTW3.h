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

namespace Physica::Core::Internal {
    template<class ScalarType> class FFTImpl {
        using RealType = typename ScalarType::ScalarType;
        using ComplexType = ComplexScalar<RealType>;
        static constexpr bool isComplex = ScalarType::isComplex;
    private:
        Vector<ComplexType> data;
        ScalarType distance;
    public:
        FFTImpl(const Vector<ScalarType>& data_, const ScalarType& distance_);
        FFTImpl(const FFTImpl& fft);
        FFTImpl(FFTImpl&& fft) noexcept;
        ~FFTImpl() = default;
        /* Operators */
        FFTImpl& operator=(const FFTImpl&) = delete;
        FFTImpl& operator=(FFTImpl&&) noexcept = delete;
        /* Operations */
        void transform();
        void invTransform();
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return data.getLength(); }
        [[nodiscard]] const ScalarType& getDistance() const noexcept { return distance; }
        [[nodiscard]] ComplexType getComponent(ssize_t index) const;
        [[nodiscard]] Vector<ComplexType> getComponents() const;
        /* Helpers */
        void swap(FFTImpl& fft);
    private:
        void transformImpl(const RealType& phase);
    };

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(const Vector<ScalarType>& data_, const ScalarType& distance_)
            : data(data_), distance(distance_) {
        assert(data.getLength() % 2 == 0);
    }

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(const FFTImpl& fft) : data(fft.data), distance(fft.distance) {}

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(FFTImpl&& fft) noexcept : data(std::move(fft.data)), distance(std::move(fft.distance)) {}

    template<class ScalarType>
    inline void FFTImpl<ScalarType>::transform() {
        transformImpl(RealType(-2 * M_PI / data.getLength()));
        data *= distance;
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType>::invTransform() {
        transformImpl(RealType(2 * M_PI / data.getLength()));
        data /= distance * data.getLength();
    }

    template<class ScalarType>
    void FFTImpl<ScalarType>::swap(FFTImpl& fft) {
        swap(data, fft.data);
        swap(distance, fft.distance);
    }

    template<class ScalarType>
    typename FFTImpl<ScalarType>::ComplexType FFTImpl<ScalarType>::getComponent(ssize_t index) const {
        const size_t length = data.getLength();
        assert(length <= SSIZE_MAX);
        assert(index <= static_cast<ssize_t>(length) / 2);
        assert(-static_cast<ssize_t>(length) / 2 <= index);
        if (index < 0)
            index += static_cast<ssize_t>(length);
        return data[index];
    }

    template<class ScalarType>
    Vector<typename FFTImpl<ScalarType>::ComplexType> FFTImpl<ScalarType>::getComponents() const {
        Vector<ComplexType> result = Vector<ComplexType>(data.getLength());
        const ssize_t half_size = static_cast<ssize_t>(data.getLength() / 2);
        for (ssize_t i = -half_size; i < half_size; ++i)
            result[i + half_size] = getComponent(i);
        return result;
    }

    template<class ScalarType>
    void FFTImpl<ScalarType>::transformImpl(const RealType& phase) {
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
