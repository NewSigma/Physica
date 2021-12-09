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

#include <fftw3.h>

namespace Physica::Core::Internal {
    template<class ScalarType> class FFTImpl {
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
        static constexpr bool isComplex = ScalarType::isComplex;
    private:
        fftw_plan forward_plan;
        fftw_plan backward_plan;
        union {
            double* real_buffer;
            fftw_complex* buffer;
        };
        int size;
        RealType deltaT;
    public:
        FFTImpl(const Vector<ScalarType>& data_, const RealType& deltaT_);
        FFTImpl(const FFTImpl& fft);
        FFTImpl(FFTImpl&& fft) noexcept;
        ~FFTImpl();
        /* Operators */
        FFTImpl& operator=(const FFTImpl&) = delete;
        FFTImpl& operator=(FFTImpl&&) noexcept = delete;
        /* Operations */
        void transform();
        void invTransform();
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return size; }
        [[nodiscard]] const RealType& getDeltaT() const noexcept { return deltaT; }
        [[nodiscard]] ComplexType getComponent(ssize_t index) const;
        [[nodiscard]] Vector<ComplexType> getComponents() const;
        /* Helpers */
        void swap(FFTImpl& fft);
    };

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(const Vector<ScalarType>& data_, const RealType& deltaT_)
            : size(static_cast<int>(data_.getLength()))
            , deltaT(deltaT_) {
        assert(data_.getLength() <= INT_MAX);
        assert(size % 2 == 0);

        if constexpr (isComplex) {
            buffer = reinterpret_cast<fftw_complex*>(fftw_malloc(data_.getLength() * sizeof(fftw_complex)));
            forward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            for (int i = 0; i < size; ++i) {
                const auto& complex = data_[i];
                buffer[i][0] = double(complex.getReal());
                buffer[i][1] = double(complex.getImag());
            }
        }
        else {
            buffer = reinterpret_cast<fftw_complex*>(fftw_malloc((size / 2 + 1) * sizeof(fftw_complex)));
            forward_plan = fftw_plan_dft_r2c_1d(size, real_buffer, buffer, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_c2r_1d(size, buffer, real_buffer, FFTW_ESTIMATE);
            for (int i = 0; i < size; ++i)
                real_buffer[i] = double(data_[i]);
        }
    }

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(const FFTImpl& fft)
            : buffer(fftw_malloc(fft.size * sizeof(fftw_complex)))
            , size(fft.size)
            , deltaT(fft.deltaT) {
        forward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    template<class ScalarType>
    FFTImpl<ScalarType>::FFTImpl(FFTImpl&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , size(fft.size)
            , deltaT(std::move(fft.deltaT)) {
        fft.plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType>
    FFTImpl<ScalarType>::~FFTImpl() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(buffer);
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType>::transform() {
        fftw_execute(forward_plan);
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType>::invTransform() {
        fftw_execute(backward_plan);
        const double factor = 1.0 / size;
        for (int i = 0; i < size; ++i) {
            buffer[i][0] *= factor;
            if constexpr (isComplex)
                buffer[i][1] *= factor;
        }
    }

    template<class ScalarType>
    void FFTImpl<ScalarType>::swap(FFTImpl& fft) {
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        std::swap(size, fft.size);
        deltaT.swap(fft.deltaT);
    }

    template<class ScalarType>
    typename FFTImpl<ScalarType>::ComplexType FFTImpl<ScalarType>::getComponent(ssize_t index) const {
        assert(index <= size / 2);
        assert(-size / 2 <= index);
        if constexpr (isComplex) {
            if (index < 0)
                index += size;
        }
        else
            index = std::abs(index);
        return ComplexType(RealType(buffer[index][0]), RealType(buffer[index][1])) * deltaT;
    }

    template<class ScalarType>
    Vector<typename FFTImpl<ScalarType>::ComplexType> FFTImpl<ScalarType>::getComponents() const {
        if constexpr (isComplex) {
            const int result_size = size + 1;
            Vector<ComplexType> result = Vector<ComplexType>(result_size);
            const ssize_t half_size = static_cast<ssize_t>(size / 2);
            for (ssize_t i = -half_size; i <= half_size; ++i)
                result[i + half_size] = getComponent(i);
            return result;
        }
        else {
            const int result_size = size / 2 + 1;
            Vector<ComplexType> result = Vector<ComplexType>(result_size);
            for (ssize_t i = 0; i < result_size; ++i)
                result[i] = getComponent(i);
            return result;
        }
    }
}
