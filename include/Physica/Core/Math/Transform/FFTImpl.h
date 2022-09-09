/*
 * Copyright 2020-2022 WeiBo He.
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
    template<class ScalarType>
    FFT<ScalarType, 1>::FFT() : forward_plan(nullptr), backward_plan(nullptr), buffer(nullptr), rSpaceSize(0), rSpaceDelta() {}

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(size_t rSpaceSize_, const RealType& rSpaceDelta_)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , rSpaceSize(static_cast<int>(rSpaceSize_))
            , rSpaceDelta(rSpaceDelta_) {
        assert(rSpaceSize_ <= static_cast<size_t>(std::numeric_limits<int>::max()));

        buffer = reinterpret_cast<fftw_complex*>(fftw_malloc(getKSpaceSize() * sizeof(fftw_complex)));
        initializePlan();
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(const Vector<ScalarType>& data_, const RealType& rSpaceDelta_)
            : FFT(data_.getLength(), rSpaceDelta_) {
        transform(data_);
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(const FFT& fft)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(fftw_malloc(fft.rSpaceSize * sizeof(fftw_complex)))
            , rSpaceSize(fft.rSpaceSize)
            , rSpaceDelta(fft.rSpaceDelta) {
        initializePlan();
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(fft.rSpaceSize)
            , rSpaceDelta(std::move(fft.rSpaceDelta)) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::~FFT() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(buffer);
    }

    template<class ScalarType>
    FFT<ScalarType, 1>& FFT<ScalarType, 1>::operator=(FFT<ScalarType, 1> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::transform(const Vector<ScalarType>& data) {
        assert(data.getLength() == static_cast<size_t>(rSpaceSize));
        if constexpr (isComplex) {
            for (int i = 0; i < rSpaceSize; ++i) {
                const auto& complex = data[i];
                buffer[i][0] = double(complex.getReal());
                buffer[i][1] = double(complex.getImag());
            }
        }
        else {
            for (int i = 0; i < rSpaceSize; ++i)
                real_buffer[i] = double(data[i]);
        }
        transform();
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::invTransform(const Vector<ComplexType>& data) {
        const size_t kSpaceSize = getKSpaceSize();
        assert(data.getLength() == kSpaceSize);
        for (size_t i = 0; i < kSpaceSize; ++i) {
            const auto& complex = data[i];
            buffer[i][0] = double(complex.getReal());
            buffer[i][1] = double(complex.getImag());
        }
        invTransform();
    }

    template<class ScalarType>
    Vector<ScalarType> FFT<ScalarType, 1>::getRSpace() const {
        Vector<ScalarType> result = Vector<ScalarType>(rSpaceSize);
        for (ssize_t i = 0; i < rSpaceSize; ++i) {
            if constexpr (isComplex)
                result[i] = ComplexType(buffer[i][0], buffer[i][1]);
            else
                result[i] = RealType(real_buffer[i]);
        }
        return result;
    }

    template<class ScalarType>
    Vector<typename FFT<ScalarType, 1>::ComplexType> FFT<ScalarType, 1>::getKSpace() const {
        const int kSpaceSize = getKSpaceSize();
        Vector<ComplexType> result = Vector<ComplexType>(kSpaceSize);
        for (int i = 0; i < kSpaceSize; ++i)
            result[i] = ComplexType(RealType(buffer[i][0]), RealType(buffer[i][1]));
        return result;
    }

    template<class ScalarType>
    void FFT<ScalarType, 1>::swap(FFT& fft) noexcept {
        using std::swap;
        swap(forward_plan, fft.forward_plan);
        swap(backward_plan, fft.backward_plan);
        swap(buffer, fft.buffer);
        swap(rSpaceSize, fft.rSpaceSize);
        swap(rSpaceDelta, fft.rSpaceDelta);
    }

    template<class ScalarType>
    int FFT<ScalarType, 1>::rSizeToKSize(int size_data) noexcept {
        if constexpr (isComplex)
            return size_data;
        else
            return size_data / 2 + 1;
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::initializePlan() {
        assert(forward_plan == nullptr);
        assert(backward_plan == nullptr);
        if constexpr (isComplex) {
            forward_plan = fftw_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        else {
            forward_plan = fftw_plan_dft_r2c_1d(rSpaceSize, real_buffer, buffer, FFTW_ESTIMATE);
            backward_plan = fftw_plan_dft_c2r_1d(rSpaceSize, buffer, real_buffer, FFTW_ESTIMATE);
        }
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::transform() {
        fftw_execute(forward_plan);
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::invTransform() {
        fftw_execute(backward_plan);
        const double factor = 1.0 / rSpaceSize;
        for (int i = 0; i < rSpaceSize; ++i) {
            if constexpr (isComplex) {
                buffer[i][0] *= factor;
                buffer[i][1] *= factor;
            }
            else
                real_buffer[i] *= factor;
        }
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT() : forward_plan(nullptr), backward_plan(nullptr), buffer(nullptr), rSpaceSize(Dim, 0), kSpaceSize(Dim, 0), rSpaceDelta() {}

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_)
            : rSpaceSize(rSpaceSize_.getLength()), rSpaceDelta(std::move(rSpaceDelta_)) {
        assert(checkSize(rSpaceSize_));
        for (size_t i = 0; i < rSpaceSize_.getLength(); ++i)
            rSpaceSize[i] = static_cast<int>(rSpaceSize_[i]);
        kSpaceSize = makeKSpaceSize();

        buffer = reinterpret_cast<fftw_complex*>(fftw_malloc(totalKSpaceSize(0) * sizeof(fftw_complex)));
        forward_plan = forwardPlan();
        backward_plan = backwardPlan();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(
            const Vector<ScalarType>& data,
            const Utils::Array<size_t, Dim>& rSpaceSize_,
            Utils::Array<RealType, Dim> rSpaceDelta_) : FFT(rSpaceSize_, std::move(rSpaceDelta_)) {
        transform(data);
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const FFT& fft)
            : buffer(fftw_malloc(fft.totalKSpaceSize(0) * sizeof(fftw_complex)))
            , rSpaceSize(fft.rSpaceSize)
            , kSpaceSize(fft.kSpaceSize)
            , rSpaceDelta(fft.rSpaceDelta) {
        forward_plan = forwardPlan();
        backward_plan = backwardPlan();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(std::move(fft.rSpaceSize))
            , kSpaceSize(std::move(fft.kSpaceSize))
            , rSpaceDelta(std::move(fft.rSpaceDelta)) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::~FFT() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(buffer);
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>& FFT<ScalarType, Dim>::operator=(FFT<ScalarType, Dim> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::swap(FFT& fft) noexcept {
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        rSpaceSize.swap(fft.rSpaceSize);
        kSpaceSize.swap(fft.kSpaceSize);
        rSpaceDelta.swap(fft.rSpaceDelta);
    }

    template<class ScalarType, size_t Dim>
    bool FFT<ScalarType, Dim>::checkSize(const Utils::Array<size_t, Dim>& rSpaceSize) {
        for (size_t elem : rSpaceSize)
            if (elem > static_cast<size_t>(std::numeric_limits<int>::max()))
                return false;
        return true;
    }

    template<class ScalarType, size_t Dim>
    Utils::Array<int, Dim> FFT<ScalarType, Dim>::makeKSpaceSize() const {
        Utils::Array<int, Dim> result(getDimen());
        size_t i = 0;
        for (; i < getDimen() - 1; ++i)
            result[i] = getRSpaceSize()[i];
        result[i] = FFT<ScalarType>::rSizeToKSize(getRSpaceSize()[i]);
        return result;
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::transform(const Vector<ScalarType>& data) {
        assert(data.getLength() == totalRSpaceSize(0));
        if constexpr (isComplex) {
            for (size_t i = 0; i < data.getLength(); ++i) {
                buffer[i][0] = double(data[i].getReal());
                buffer[i][1] = double(data[i].getImag());
            }
        }
        else {
            const size_t lastDimSize = getRSpaceSize(getDimen() - 1);
            const size_t paddingCount = (lastDimSize / 2 + 1) * 2 - lastDimSize;
            size_t j = 0;
            for (size_t i = 0; i < data.getLength(); ++i, ++j) {
                if (i % lastDimSize == 0 && i != 0)
                    j += paddingCount;
                real_buffer[j] = double(data[i]);
            }
        }
        transform();
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::invTransform(const Vector<ComplexType>& data) {
        assert(data.getLength() == totalKSpaceSize(0));
        for (size_t i = 0; i < data.getLength(); ++i) {
            const auto& complex = data[i];
            buffer[i][0] = double(complex.getReal());
            buffer[i][1] = double(complex.getImag());
        }
        invTransform();
    }

    template<class ScalarType, size_t Dim>
    Vector<ScalarType> FFT<ScalarType, Dim>::getRSpace() const {
        Vector<ScalarType> result = Vector<ScalarType>(totalRSpaceSize(0));
        if constexpr (isComplex) {
            for (size_t i = 0; i < result.getLength(); ++i)
                result[i] = ComplexType(buffer[i][0], buffer[i][1]);
        }
        else {
            const size_t lastDimSize = getRSpaceSize(getDimen() - 1);
            const size_t paddingCount = (lastDimSize / 2 + 1) * 2 - lastDimSize;
            size_t j = 0;
            for (size_t i = 0; i < result.getLength(); ++i, ++j) {
                if (i % lastDimSize == 0 && i != 0)
                    j += paddingCount;
                result[i] = RealType(real_buffer[j]);
            }
        }    
        return result;
    }

    template<class ScalarType, size_t Dim>
    Vector<typename FFT<ScalarType, Dim>::ComplexType> FFT<ScalarType, Dim>::getKSpace() const {
        Vector<ComplexType> result = Vector<ComplexType>(totalKSpaceSize(0));
        for (size_t i = 0; i < result.getLength(); ++i)
            result[i] = ComplexType(RealType(buffer[i][0]), RealType(buffer[i][1]));
        return result;
    }

    template<class ScalarType, size_t Dim>
    fftw_plan FFT<ScalarType, Dim>::forwardPlan() {
        if constexpr (Dim == 2)
            if constexpr (isComplex)
                return fftw_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c_2d(rSpaceSize[0], rSpaceSize[1], real_buffer, buffer, FFTW_ESTIMATE);
        else if constexpr (Dim == 3)
            if constexpr (isComplex)
                return fftw_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], real_buffer, buffer, FFTW_ESTIMATE);
        else
            if constexpr (isComplex)
                return fftw_plan_dft(getDimen(), rSpaceSize.data(), buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_r2c(getDimen(), rSpaceSize.data(), real_buffer, buffer, FFTW_ESTIMATE);
    }
    
    template<class ScalarType, size_t Dim>
    fftw_plan FFT<ScalarType, Dim>::backwardPlan() {
        if constexpr (Dim == 2)
            if constexpr (isComplex)
                return fftw_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r_2d(rSpaceSize[0], rSpaceSize[1], buffer, real_buffer, FFTW_ESTIMATE);
        else if constexpr (Dim == 3)
            if constexpr (isComplex)
                return fftw_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, real_buffer, FFTW_ESTIMATE);
        else
            if constexpr (isComplex)
                return fftw_plan_dft(getDimen(), rSpaceSize.data(), buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            else
                return fftw_plan_dft_c2r(getDimen(), rSpaceSize.data(), buffer, real_buffer, FFTW_ESTIMATE);
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::transform() {
        fftw_execute(forward_plan);
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::invTransform() {
        fftw_execute(backward_plan);
        const size_t totalSize = totalRSpaceSize(0);
        const double factor = 1.0 / totalSize;
        if constexpr (isComplex) {
            for (size_t i = 0; i < totalSize; ++i) {
                buffer[i][0] *= factor;
                buffer[i][1] *= factor;
            }
        }
        else {
            const size_t lastDimSize = getRSpaceSize(getDimen() - 1);
            const size_t paddingCount = (lastDimSize / 2 + 1) * 2 - lastDimSize;
            size_t j = 0;
            for (size_t i = 0; i < totalSize; ++i, ++j) {
                if (i % lastDimSize == 0 && i != 0)
                    j += paddingCount;
                real_buffer[j] *= factor;
            }
        }
    }

    template<class ScalarType, size_t Dim>
    size_t FFT<ScalarType, Dim>::totalRSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDimen(); ++i)
            result *= getRSpaceSize(i);
        return result;
    }

    template<class ScalarType, size_t Dim>
    size_t FFT<ScalarType, Dim>::totalKSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDimen(); ++i)
            result *= getKSpaceSize(i);
        return result;
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const {
        for (size_t i = 0; i < getDimen(); ++i) {
            const int size_i = rSpaceSize[i];
            ssize_t index = indexes[i];
            assert(index <= size_i / 2);
            assert(-size_i / 2 <= index);
            if (index < 0)
                index += size_i;
            indexes[i] = index;
        }
    }

    template<class ScalarType, size_t Dim>
    typename FFT<ScalarType, Dim>::RealType FFT<ScalarType, Dim>::mulDeltas() const {
        RealType result = rSpaceDelta[0];
        for (size_t i = 1; i < rSpaceDelta.getLength(); ++i)
            result *= rSpaceDelta[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    size_t FFT<ScalarType, Dim>::componentsSizeFrom(size_t dim) const {
        size_t result = 1;
        for (size_t i = dim; i < getDimen(); ++i) {
            if constexpr (isComplex)
                result *= rSpaceSize[i] / 2 * 2 + 1;
            else {
                if (i == getDimen() - 1)
                    result *= rSpaceSize[i] / 2 + 1;
                else
                    result *= rSpaceSize[i] / 2 * 2 + 1;
            }
                
        }
        return result;
    }

    template<class ScalarType, size_t Dim>
    Utils::Array<ssize_t, Dim> FFT<ScalarType, Dim>::linearIndexToDim(size_t index) const {
        Utils::Array<ssize_t, Dim> result(getDimen());
        for (size_t i = 0; i < getDimen(); ++i) {
            const size_t componentsSizeFrom_i = componentsSizeFrom(i + 1);
            ssize_t dim_i = index / componentsSizeFrom_i;
            index -= componentsSizeFrom_i * dim_i;
            if constexpr (isComplex)
                result[i] = dim_i - rSpaceSize[i] / 2;
            else {
                if (i == getDimen() - 1)
                    result[i] = dim_i;
                else
                    result[i] = dim_i - rSpaceSize[i] / 2;
            }
        }
        return result;
    }
}
