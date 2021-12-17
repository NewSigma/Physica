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
    template<class ScalarType, size_t Dim> class FFTImpl;

    template<class ScalarType>
    class FFTImpl<ScalarType, 1> {
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
        FFTImpl(size_t size_, const RealType& deltaT_);
        FFTImpl(const Vector<ScalarType>& data_, const RealType& deltaT_);
        FFTImpl(const FFTImpl& fft);
        FFTImpl(FFTImpl&& fft) noexcept;
        ~FFTImpl();
        /* Operators */
        FFTImpl& operator=(const FFTImpl&) = delete;
        FFTImpl& operator=(FFTImpl&&) noexcept = delete;
        /* Operations */
        void transform(const Vector<ScalarType>& data);
        void invTransform(const Vector<ComplexType>& data);
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return size; }
        [[nodiscard]] const RealType& getDeltaT() const noexcept { return deltaT; }
        [[nodiscard]] ComplexType getComponent(ssize_t index) const;
        [[nodiscard]] Vector<ComplexType> getComponents() const;
        /* Helpers */
        void swap(FFTImpl& fft);
    private:
        void transform();
        void invTransform();
    };

    template<class ScalarType>
    FFTImpl<ScalarType, 1>::FFTImpl(size_t size_, const RealType& deltaT_)
            : size(static_cast<int>(size_))
            , deltaT(deltaT_) {
        assert(size_ <= INT_MAX);

        if constexpr (isComplex) {
            buffer = reinterpret_cast<fftw_complex*>(fftw_malloc(size * sizeof(fftw_complex)));
            forward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        else {
            buffer = reinterpret_cast<fftw_complex*>(fftw_malloc((size / 2 + 1) * sizeof(fftw_complex)));
            forward_plan = fftw_plan_dft_r2c_1d(size, real_buffer, buffer, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_c2r_1d(size, buffer, real_buffer, FFTW_ESTIMATE);
        }
    }

    template<class ScalarType>
    FFTImpl<ScalarType, 1>::FFTImpl(const Vector<ScalarType>& data_, const RealType& deltaT_)
            : FFTImpl(data_.getLength(), deltaT_) {
        transform(data_);
    }

    template<class ScalarType>
    FFTImpl<ScalarType, 1>::FFTImpl(const FFTImpl& fft)
            : buffer(fftw_malloc(fft.size * sizeof(fftw_complex)))
            , size(fft.size)
            , deltaT(fft.deltaT) {
        forward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_1d(size, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    template<class ScalarType>
    FFTImpl<ScalarType, 1>::FFTImpl(FFTImpl&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , size(fft.size)
            , deltaT(std::move(fft.deltaT)) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType>
    FFTImpl<ScalarType, 1>::~FFTImpl() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(buffer);
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType, 1>::transform(const Vector<ScalarType>& data) {
        assert(data.getLength() == static_cast<size_t>(size));
        if constexpr (isComplex) {
            for (int i = 0; i < size; ++i) {
                const auto& complex = data[i];
                buffer[i][0] = double(complex.getReal());
                buffer[i][1] = double(complex.getImag());
            }
        }
        else {
            for (int i = 0; i < size; ++i)
                real_buffer[i] = double(data[i]);
        }
        transform();
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType, 1>::invTransform(const Vector<ComplexType>& data) {
        [[maybe_unused]] const size_t expectedSize = isComplex ? size : size / 2 + 1;
        assert(data.getLength() == expectedSize);
        for (int i = 0; i < expectedSize; ++i) {
            const auto& complex = data[i];
            buffer[i][0] = double(complex.getReal());
            buffer[i][1] = double(complex.getImag());
        }
        invTransform();
    }

    template<class ScalarType>
    typename FFTImpl<ScalarType, 1>::ComplexType FFTImpl<ScalarType, 1>::getComponent(ssize_t index) const {
        assert(index <= size / 2);
        assert(-size / 2 <= index);
        if constexpr (isComplex) {
            if (index < 0)
                index += size;
        }
        else
            assert(index >= 0);
        return ComplexType(RealType(buffer[index][0]), RealType(buffer[index][1])) * deltaT;
    }

    template<class ScalarType>
    Vector<typename FFTImpl<ScalarType, 1>::ComplexType> FFTImpl<ScalarType, 1>::getComponents() const {
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

    template<class ScalarType>
    void FFTImpl<ScalarType, 1>::swap(FFTImpl& fft) {
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        std::swap(size, fft.size);
        std::swap(deltaT, fft.deltaT);
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType, 1>::transform() {
        fftw_execute(forward_plan);
    }

    template<class ScalarType>
    inline void FFTImpl<ScalarType, 1>::invTransform() {
        fftw_execute(backward_plan);
        const double factor = 1.0 / size;
        for (int i = 0; i < size; ++i) {
            if constexpr (isComplex) {
                buffer[i][0] *= factor;
                buffer[i][1] *= factor;
            }
            else
                real_buffer[i] *= factor;
        }
    }

    template<class ScalarType, size_t Dim>
    class FFTImpl {
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
        Utils::Array<int, Dim> size;
        Utils::Array<RealType, Dim> deltaTs;
    public:
        FFTImpl(const Utils::Array<size_t, Dim>& size_,
                Utils::Array<RealType, Dim> deltaTs_);
        FFTImpl(const Vector<ScalarType>& data,
                const Utils::Array<size_t, Dim>& size_,
                Utils::Array<RealType, Dim> deltaTs_);
        FFTImpl(const FFTImpl& fft);
        FFTImpl(FFTImpl&& fft) noexcept;
        ~FFTImpl();
        /* Operators */
        FFTImpl& operator=(const FFTImpl&) = delete;
        FFTImpl& operator=(FFTImpl&&) noexcept = delete;
        /* Operations */
        void transform(const Vector<ScalarType>& data);
        void invTransform(const Vector<ComplexType>& data);
        /* Getters */
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return size[dim]; }
        [[nodiscard]] const RealType& getDeltaT(size_t dim) const noexcept { return deltaTs[dim]; }
        [[nodiscard]] size_t getDimen() const noexcept { return size.getLength(); }
        [[nodiscard]] ComplexType getComponent(Utils::Array<ssize_t, Dim> indexes) const;
        [[nodiscard]] Vector<ComplexType> getComponents() const;
        /* Helpers */
        void swap(FFTImpl& fft);
    private:
        void transform();
        void invTransform();
        [[nodiscard]] fftw_plan forwardPlan();
        [[nodiscard]] fftw_plan backwardPlan();
        [[nodiscard]] size_t outputSizeFrom(size_t dim) const;
        void normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const;
        [[nodiscard]] RealType mulDeltaTs() const;
        [[nodiscard]] size_t componentsSizeFrom(size_t dim) const;
        [[nodiscard]] Utils::Array<ssize_t, Dim> linearIndexToDim(size_t index) const;
    };

    template<class ScalarType, size_t Dim>
    FFTImpl<ScalarType, Dim>::FFTImpl(
            const Utils::Array<size_t, Dim>& size_,
            Utils::Array<RealType, Dim> deltaTs_) : size(size_.getLength()), deltaTs(std::move(deltaTs_)) {
        for (size_t i = 0; i < getDimen(); ++i)
            size[i] = static_cast<int>(size_[i]);

        buffer = reinterpret_cast<fftw_complex*>(fftw_malloc(outputSizeFrom(0) * sizeof(fftw_complex)));
        forward_plan = forwardPlan();
        backward_plan = backwardPlan();
    }

    template<class ScalarType, size_t Dim>
    FFTImpl<ScalarType, Dim>::FFTImpl(
            const Vector<ScalarType>& data,
            const Utils::Array<size_t, Dim>& size_,
            Utils::Array<RealType, Dim> deltaTs_) : FFTImpl(size_, std::move(deltaTs_)) {
        transform(data);
    }

    template<class ScalarType, size_t Dim>
    FFTImpl<ScalarType, Dim>::FFTImpl(const FFTImpl& fft)
            : buffer(fftw_malloc(fft.outputSizeFrom(0) * sizeof(fftw_complex)))
            , size(fft.size)
            , deltaTs(fft.deltaTs) {
        forward_plan = forwardPlan();
        backward_plan = backwardPlan();
    }

    template<class ScalarType, size_t Dim>
    FFTImpl<ScalarType, Dim>::FFTImpl(FFTImpl&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , size(std::move(fft.size))
            , deltaTs(std::move(fft.deltaTs)) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType, size_t Dim>
    FFTImpl<ScalarType, Dim>::~FFTImpl() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(buffer);
    }

    template<class ScalarType, size_t Dim>
    inline void FFTImpl<ScalarType, Dim>::transform(const Vector<ScalarType>& data) {
        if constexpr (isComplex) {
            for (size_t i = 0; i < data.getLength(); ++i) {
                buffer[i][0] = double(data[i].getReal());
                buffer[i][1] = double(data[i].getImag());
            }
        }
        else {
            const size_t lastDimSize = size[getDimen() - 1];
            const size_t paddingCount = (lastDimSize / 2 + 1) * 2 - lastDimSize;
            size_t j = 0;
            for (size_t i = 0; i < data.getLength(); ++i, ++j) {
                real_buffer[j] = double(data[i]);
                if (i % lastDimSize == 0)
                    j += paddingCount;
            }
        }
        transform();
    }

    template<class ScalarType, size_t Dim>
    inline void FFTImpl<ScalarType, Dim>::invTransform(const Vector<ComplexType>& data) {
        [[maybe_unused]] const size_t expectedSize = outputSizeFrom(0);
        assert(data.getLength() == expectedSize);
        for (int i = 0; i < expectedSize; ++i) {
            const auto& complex = data[i];
            buffer[i][0] = double(complex.getReal());
            buffer[i][1] = double(complex.getImag());
        }
        invTransform();
    }

    template<class ScalarType, size_t Dim>
    void FFTImpl<ScalarType, Dim>::swap(FFTImpl& fft) {
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        std::swap(size, fft.size);
        std::swap(deltaTs, fft.deltaTs);
    }

    template<class ScalarType, size_t Dim>
    typename FFTImpl<ScalarType, Dim>::ComplexType FFTImpl<ScalarType, Dim>::getComponent(Utils::Array<ssize_t, Dim> indexes) const {
        normalizeIndexes(indexes);
        size_t index = 0;
        for (size_t i = 0; i < indexes.getLength(); ++i)
            index += indexes[i] * outputSizeFrom(i + 1);
        return ComplexType(RealType(buffer[index][0]), RealType(buffer[index][1])) * mulDeltaTs();
    }

    template<class ScalarType, size_t Dim>
    Vector<typename FFTImpl<ScalarType, Dim>::ComplexType> FFTImpl<ScalarType, Dim>::getComponents() const {
        Vector<ComplexType> result = Vector<ComplexType>(componentsSizeFrom(0));
        for (size_t i = 0; i < result.getLength(); ++i)
            result[i] = getComponent(linearIndexToDim(i));
        return result;
    }

    template<class ScalarType, size_t Dim>
    inline void FFTImpl<ScalarType, Dim>::transform() {
        fftw_execute(forward_plan);
    }

    template<class ScalarType, size_t Dim>
    inline void FFTImpl<ScalarType, Dim>::invTransform() {
        fftw_execute(backward_plan);
        const size_t totalSize = outputSizeFrom(0);
        const double factor = 1.0 / totalSize;
        for (int i = 0; i < totalSize; ++i) {
            if constexpr (isComplex) {
                buffer[i][0] *= factor;
                buffer[i][1] *= factor;
            }
            else
                real_buffer[i] *= factor;
        }
    }

    template<class ScalarType, size_t Dim>
    fftw_plan FFTImpl<ScalarType, Dim>::forwardPlan() {
        if constexpr (Dim == 2)
            if constexpr (isComplex)
                return fftw_plan_dft_2d(size[0], size[1], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c_2d(size[0], size[1], real_buffer, buffer, FFTW_ESTIMATE);
        else if constexpr (Dim == 3)
            if constexpr (isComplex)
                return fftw_plan_dft_3d(size[0], size[1], size[2], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c_3d(size[0], size[1], size[2], real_buffer, buffer, FFTW_ESTIMATE);
        else
            if constexpr (isComplex)
                return fftw_plan_dft(getDimen(), size.data(), buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c(getDimen(), size.data(), real_buffer, buffer, FFTW_ESTIMATE);
    }
    
    template<class ScalarType, size_t Dim>
    fftw_plan FFTImpl<ScalarType, Dim>::backwardPlan() {
        if constexpr (Dim == 2)
            if constexpr (isComplex)
                return fftw_plan_dft_2d(size[0], size[1], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r_2d(size[0], size[1], buffer, real_buffer, FFTW_ESTIMATE);
        else if constexpr (Dim == 3)
            if constexpr (isComplex)
                return fftw_plan_dft_3d(size[0], size[1], size[2], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r_3d(size[0], size[1], size[2], buffer, real_buffer, FFTW_ESTIMATE);
        else
            if constexpr (isComplex)
                return fftw_plan_dft(getDimen(), size.data(), buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r(getDimen(), size.data(), buffer, real_buffer, FFTW_ESTIMATE);
    }

    template<class ScalarType, size_t Dim>
    size_t FFTImpl<ScalarType, Dim>::outputSizeFrom(size_t dim) const {
        size_t result = 1;
        for (size_t i = dim; i < getDimen(); ++i) {
            if constexpr (isComplex)
                result *= size[i];
            else {
                if (i == getDimen() - 1)
                    result *= size[i] / 2 + 1;
                else
                    result *= size[i];
            }
        }
        return result;
    }

    template<class ScalarType, size_t Dim>
    void FFTImpl<ScalarType, Dim>::normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const {
        for (size_t i = 0; i < getDimen(); ++i) {
            const int size_i = size[i];
            ssize_t index = indexes[i];
            assert(index <= size_i / 2);
            assert(-size_i / 2 <= index);
            if (index < 0)
                index += size_i;
            indexes[i] = index;
        }
    }

    template<class ScalarType, size_t Dim>
    typename FFTImpl<ScalarType, Dim>::RealType FFTImpl<ScalarType, Dim>::mulDeltaTs() const {
        RealType result = deltaTs[0];
        for (size_t i = 1; i < deltaTs.getLength(); ++i)
            result *= deltaTs[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    size_t FFTImpl<ScalarType, Dim>::componentsSizeFrom(size_t dim) const {
        size_t result = 1;
        for (size_t i = dim; i < getDimen(); ++i) {
            if constexpr (isComplex)
                result *= size[i] / 2 * 2 + 1;
            else {
                if (i == getDimen() - 1)
                    result *= size[i] / 2 + 1;
                else
                    result *= size[i] / 2 * 2 + 1;
            }
                
        }
        return result;
    }

    template<class ScalarType, size_t Dim>
    Utils::Array<ssize_t, Dim> FFTImpl<ScalarType, Dim>::linearIndexToDim(size_t index) const {
        Utils::Array<ssize_t, Dim> result(getDimen());
        for (size_t i = 0; i < getDimen(); ++i) {
            const size_t componentsSizeFrom_i = componentsSizeFrom(i + 1);
            ssize_t dim_i = index / componentsSizeFrom_i;
            index -= componentsSizeFrom_i * dim_i;
            if constexpr (isComplex)
                result[i] = dim_i - size[i] / 2;
            else {
                if (i == getDimen() - 1)
                    result[i] = dim_i;
                else
                    result[i] = dim_i - size[i] / 2;
            }
        }
        return result;
    }
}
