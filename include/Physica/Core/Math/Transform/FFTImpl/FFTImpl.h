/*
 * Copyright 2020-2023 WeiBo He.
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

#include <mutex>

namespace Physica::Core {
    namespace Internal {
        class ThreadGuardFFTW final {
        public:
            std::mutex globalMutex;
        public:
            ThreadGuardFFTW(const ThreadGuardFFTW&) = delete;
            ThreadGuardFFTW(ThreadGuardFFTW&&) noexcept = delete;
            ~ThreadGuardFFTW() = default;
            /* Operators */
            ThreadGuardFFTW& operator=(const ThreadGuardFFTW&) = delete;
            ThreadGuardFFTW& operator=(ThreadGuardFFTW&&) noexcept = delete;
            /* Static members */
            [[nodiscard]] static ThreadGuardFFTW& getInstance();
        private:
            ThreadGuardFFTW() = default;
        };
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT()
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(nullptr)
            , rSpaceSize(0)
            , rSpaceDelta()
            , planFlag(PlanFlag::Measure) {}

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(size_t rSpaceSize_, RealType rSpaceDelta_)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , rSpaceSize(static_cast<int>(rSpaceSize_))
            , rSpaceDelta(rSpaceDelta_)
            , planFlag(PlanFlag::Measure) {
        assert(rSpaceSize_ <= static_cast<size_t>(std::numeric_limits<int>::max()));
        buffer = reinterpret_cast<ComplexTypeFFTW*>(fftw_malloc(getKSpaceSize() * sizeof(ComplexTypeFFTW)));
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(size_t rSpaceSize_, RealType rSpaceDelta_, PlanFlag planFlag_)
            : FFT(rSpaceSize_, rSpaceDelta_) {
        planFlag = planFlag_;
        initializePlan();
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(const Vector<ScalarType>& data_, RealType rSpaceDelta_, PlanFlag planFlag)
            : FFT(data_.getLength(), rSpaceDelta_, planFlag) {
        transform(data_);
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(const FFT& fft)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(reinterpret_cast<ComplexTypeFFTW*>(fftw_malloc(fft.rSpaceSize * sizeof(ComplexTypeFFTW))))
            , rSpaceSize(fft.rSpaceSize)
            , rSpaceDelta(fft.rSpaceDelta)
            , planFlag(fft.planFlag) {
        initializePlan();
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(fft.rSpaceSize)
            , rSpaceDelta(std::move(fft.rSpaceDelta))
            , planFlag(fft.planFlag) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType>
    FFT<ScalarType, 1>::~FFT() {
        Internal::ThreadGuardFFTW::getInstance().globalMutex.lock();
        if constexpr (isSinglePrec) {
            fftwf_destroy_plan(forward_plan);
            fftwf_destroy_plan(backward_plan);
        }
        else {
            fftw_destroy_plan(forward_plan);
            fftw_destroy_plan(backward_plan);
        }
        fftw_free(buffer);
        Internal::ThreadGuardFFTW::getInstance().globalMutex.unlock();
    }

    template<class ScalarType>
    FFT<ScalarType, 1>& FFT<ScalarType, 1>::operator=(FFT<ScalarType, 1> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType>
    void FFT<ScalarType, 1>::swap(FFT& fft) noexcept {
        assert(this != &fft && "[Error]: Self swap is likely a bug");
        using std::swap;
        swap(forward_plan, fft.forward_plan);
        swap(backward_plan, fft.backward_plan);
        swap(buffer, fft.buffer);
        swap(rSpaceSize, fft.rSpaceSize);
        swap(rSpaceDelta, fft.rSpaceDelta);
        swap(planFlag, fft.planFlag);
    }

    template<class ScalarType>
    inline FFT<ScalarType, 1> FFT<ScalarType, 1>::makeEmptyFFT(size_t rSpaceSize, RealType rSpaceDelta) {
        return FFT<ScalarType, 1>(rSpaceSize, rSpaceDelta);
    }

    template<class ScalarType>
    template<class IndexType>
    __host__ __device__ constexpr inline IndexType FFT<ScalarType, 1>::rSizeToKSize(IndexType size_data) noexcept {
        if constexpr (isComplex)
            return size_data;
        else
            return size_data / 2 + 1;
    }

    template<class ScalarType>
    void FFT<ScalarType, 1>::transform(const This& planProvider, This& bufferProvider) {
        const auto forward_plan = planProvider.forward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(forward_plan != nullptr);
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isSinglePrec) {
            if constexpr (isComplex)
                fftwf_execute_dft(forward_plan, buffer, buffer);
            else
                fftwf_execute_dft_r2c(forward_plan, reinterpret_cast<float*>(buffer), buffer);
        }
        else {
            if constexpr (isComplex)
                fftw_execute_dft(forward_plan, buffer, buffer);
            else
                fftw_execute_dft_r2c(forward_plan, reinterpret_cast<double*>(buffer), buffer);
        }
    }

    template<class ScalarType>
    void FFT<ScalarType, 1>::rawInvTransform(const This& planProvider, This& bufferProvider) {
        const auto backward_plan = planProvider.backward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(backward_plan != nullptr);
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isSinglePrec) {
            if constexpr (isComplex)
                fftwf_execute_dft(backward_plan, buffer, buffer);
            else
                fftwf_execute_dft_c2r(backward_plan, buffer, reinterpret_cast<float*>(buffer));
        }
        else {
            if constexpr (isComplex)
                fftw_execute_dft(backward_plan, buffer, buffer);
            else
                fftw_execute_dft_c2r(backward_plan, buffer, reinterpret_cast<double*>(buffer));
        }
    }

    template<class ScalarType>
    void FFT<ScalarType, 1>::invTransform(const This& planProvider, This& bufferProvider) {
        rawInvTransform(planProvider, bufferProvider);
        const ScalarType factor = RealType(1.0 / planProvider.getRSpaceSize());
        bufferProvider.getRSpace() *= factor;
    }

    template<class ScalarType>
    inline void FFT<ScalarType, 1>::initializePlan() {
        assert(forward_plan == nullptr);
        assert(backward_plan == nullptr);
        Internal::ThreadGuardFFTW::getInstance().globalMutex.lock();
        const int flag = static_cast<int>(planFlag);
        if constexpr (isComplex) {
            if constexpr (isSinglePrec) {
                forward_plan = fftwf_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_FORWARD, flag);
                backward_plan = fftwf_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_BACKWARD, flag);
            }
            else {
                forward_plan = fftw_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_FORWARD, flag);
                backward_plan = fftw_plan_dft_1d(rSpaceSize, buffer, buffer, FFTW_BACKWARD, flag);
            }
        }
        else {
            if constexpr (isSinglePrec) {
                forward_plan = fftwf_plan_dft_r2c_1d(rSpaceSize, (float*)buffer, buffer, flag);
                backward_plan = fftwf_plan_dft_c2r_1d(rSpaceSize, buffer, (float*)buffer, flag);
            }
            else {
                forward_plan = fftw_plan_dft_r2c_1d(rSpaceSize, (double*)buffer, buffer, flag);
                backward_plan = fftw_plan_dft_c2r_1d(rSpaceSize, buffer, (double*)buffer, flag);
            }
        }
        Internal::ThreadGuardFFTW::getInstance().globalMutex.unlock();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT()
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(nullptr)
            , rSpaceSize(Dim, 0)
            , kSpaceSize(Dim, 0)
            , rSpaceDelta()
            , planFlag(PlanFlag::Measure) {}

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , rSpaceSize(rSpaceSize_.getLength())
            , rSpaceDelta(std::move(rSpaceDelta_))
            , planFlag(PlanFlag::Measure) {
        assert(checkSize(rSpaceSize_));
        for (size_t i = 0; i < rSpaceSize_.getLength(); ++i)
            rSpaceSize[i] = static_cast<int>(rSpaceSize_[i]);
        kSpaceSize = rSizeToKSize(rSpaceSize);

        buffer = reinterpret_cast<ComplexTypeFFTW*>(fftw_malloc(sumKSpaceSize(0) * sizeof(ComplexTypeFFTW)));
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_, PlanFlag planFlag_)
            : FFT(rSpaceSize_, std::move(rSpaceDelta_)) {
        planFlag = planFlag_;
        initializePlan();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const FFT& fft)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(reinterpret_cast<ComplexTypeFFTW*>(fftw_malloc(fft.sumKSpaceSize(0) * sizeof(ComplexTypeFFTW))))
            , rSpaceSize(fft.rSpaceSize)
            , kSpaceSize(fft.kSpaceSize)
            , rSpaceDelta(fft.rSpaceDelta)
            , planFlag(fft.planFlag) {
        initializePlan();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(std::move(fft.rSpaceSize))
            , kSpaceSize(std::move(fft.kSpaceSize))
            , rSpaceDelta(std::move(fft.rSpaceDelta))
            , planFlag(fft.planFlag) {
        fft.forward_plan = nullptr;
        fft.backward_plan = nullptr;
        fft.buffer = nullptr;
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::~FFT() {
        Internal::ThreadGuardFFTW::getInstance().globalMutex.lock();
        if constexpr (isSinglePrec) {
            fftwf_destroy_plan(forward_plan);
            fftwf_destroy_plan(backward_plan);
        }
        else {
            fftw_destroy_plan(forward_plan);
            fftw_destroy_plan(backward_plan);
        }
        fftw_free(buffer);
        Internal::ThreadGuardFFTW::getInstance().globalMutex.unlock();
    }

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>& FFT<ScalarType, Dim>::operator=(FFT<ScalarType, Dim> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::swap(FFT& fft) noexcept {
        assert(this != &fft && "[Error]: Self swap is likely a bug");
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        rSpaceSize.swap(fft.rSpaceSize);
        kSpaceSize.swap(fft.kSpaceSize);
        rSpaceDelta.swap(fft.rSpaceDelta);
        std::swap(planFlag, fft.planFlag);
    }

    template<class ScalarType, size_t Dim>
    typename FFT<ScalarType, Dim>::IndexArray FFT<ScalarType, Dim>::getRSpaceSize() const noexcept{
        IndexArray result{};
        for (size_t i = 0; i < Dim; ++i)
            result[i] = rSpaceSize[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    typename FFT<ScalarType, Dim>::IndexArray FFT<ScalarType, Dim>::getKSpaceSize() const noexcept {
        IndexArray result{};
        for (size_t i = 0; i < Dim; ++i)
            result[i] = kSpaceSize[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    inline FFT<ScalarType, Dim> FFT<ScalarType, Dim>::makeEmptyFFT(
            const Utils::Array<size_t, Dim>& rSpaceSize,
            Utils::Array<RealType, Dim> rSpaceDelta) {
        return FFT(rSpaceSize, rSpaceDelta);
    }

    template<class ScalarType, size_t Dim>
    template<class IndexType>
    Utils::Array<IndexType, Dim> FFT<ScalarType, Dim>::rSizeToKSize(const Utils::Array<IndexType, Dim>& rSize) {
        Utils::Array<IndexType, Dim> result(rSize.getLength());
        size_t i = 0;
        for (; i < rSize.getLength() - 1; ++i)
            result[i] = rSize[i];
        result[i] = FFT<ScalarType>::rSizeToKSize(rSize[i]);
        return result;
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::transform(const This& planProvider, This& bufferProvider) {
        const auto forward_plan = planProvider.forward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(forward_plan != nullptr);
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isSinglePrec) {
            if constexpr (isComplex)
                fftwf_execute_dft(forward_plan, buffer, buffer);
            else
                fftwf_execute_dft_r2c(forward_plan, reinterpret_cast<float*>(buffer), buffer);
        }
        else {
            if constexpr (isComplex)
                fftw_execute_dft(forward_plan, buffer, buffer);
            else
                fftw_execute_dft_r2c(forward_plan, reinterpret_cast<double*>(buffer), buffer);
        }
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::rawInvTransform(const This& planProvider, This& bufferProvider) {
        const auto backward_plan = planProvider.backward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(backward_plan != nullptr);
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isSinglePrec) {
            if constexpr (isComplex)
                fftwf_execute_dft(backward_plan, buffer, buffer);
            else
                fftwf_execute_dft_c2r(backward_plan, buffer, reinterpret_cast<float*>(buffer));
        }
        else {
            if constexpr (isComplex)
                fftw_execute_dft(backward_plan, buffer, buffer);
            else
                fftw_execute_dft_c2r(backward_plan, buffer, reinterpret_cast<double*>(buffer));
        }
        
    }

    template<class ScalarType, size_t Dim>
    inline void FFT<ScalarType, Dim>::invTransform(const This& planProvider, This& bufferProvider) {
        rawInvTransform(planProvider, bufferProvider);
        const ScalarType factor = RealType(1.0 / planProvider.sumRSpaceSize(0));
        bufferProvider.getRSpace() *= factor;
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::initializePlan() {
        assert(forward_plan == nullptr);
        assert(backward_plan == nullptr);
        Internal::ThreadGuardFFTW::getInstance().globalMutex.lock();
        forward_plan = makeForwardPlan();
        backward_plan = makeBackwardPlan();
        Internal::ThreadGuardFFTW::getInstance().globalMutex.unlock();
    }

    template<class ScalarType, size_t Dim>
    fftw_plan FFT<ScalarType, Dim>::makeForwardPlan() {
        fftw_plan plan;
        if constexpr (Dim == 2) {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_r2c_2d(rSpaceSize[0], rSpaceSize[1], reinterpret_cast<float*>(buffer), buffer, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_r2c_2d(rSpaceSize[0], rSpaceSize[1], reinterpret_cast<double*>(buffer), buffer, FFTW_ESTIMATE);
            }
        }
        else if constexpr (Dim == 3) {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_r2c_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], reinterpret_cast<float*>(buffer), buffer, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_r2c_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], reinterpret_cast<double*>(buffer), buffer, FFTW_ESTIMATE);
            }
        }
        else {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft(getDim(), rSpaceSize.data(), buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft(getDim(), rSpaceSize.data(), buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_r2c(getDim(), rSpaceSize.data(), reinterpret_cast<float*>(buffer), buffer, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_r2c(getDim(), rSpaceSize.data(), reinterpret_cast<double*>(buffer), buffer, FFTW_ESTIMATE);
            }
        }
        return plan;
    }
    
    template<class ScalarType, size_t Dim>
    fftw_plan FFT<ScalarType, Dim>::makeBackwardPlan() {
        fftw_plan plan;
        if constexpr (Dim == 2) {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_2d(rSpaceSize[0], rSpaceSize[1], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_c2r_2d(rSpaceSize[0], rSpaceSize[1], buffer, reinterpret_cast<float*>(buffer), FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_c2r_2d(rSpaceSize[0], rSpaceSize[1], buffer, reinterpret_cast<double*>(buffer), FFTW_ESTIMATE);
            }
        }
        else if constexpr (Dim == 3) {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_c2r_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, reinterpret_cast<float*>(buffer), FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_c2r_3d(rSpaceSize[0], rSpaceSize[1], rSpaceSize[2], buffer, reinterpret_cast<double*>(buffer), FFTW_ESTIMATE);
            }
        }
        else {
            if constexpr (isComplex) {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft(getDim(), rSpaceSize.data(), buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft(getDim(), rSpaceSize.data(), buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            else {
                if constexpr (isSinglePrec)
                    plan = fftwf_plan_dft_c2r(getDim(), rSpaceSize.data(), buffer, reinterpret_cast<float*>(buffer), FFTW_ESTIMATE);
                else
                    plan = fftw_plan_dft_c2r(getDim(), rSpaceSize.data(), buffer, reinterpret_cast<double*>(buffer), FFTW_ESTIMATE);
            }
        }
        return plan;
    }

    template<class ScalarType, size_t Dim>
    size_t FFT<ScalarType, Dim>::sumRSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDim(); ++i)
            result *= getRSpaceSize()[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    size_t FFT<ScalarType, Dim>::sumKSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDim(); ++i)
            result *= getKSpaceSize()[i];
        return result;
    }

    template<class ScalarType, size_t Dim>
    void FFT<ScalarType, Dim>::normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const {
        for (size_t i = 0; i < getDim(); ++i) {
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
        for (size_t i = dim; i < getDim(); ++i) {
            if constexpr (isComplex)
                result *= rSpaceSize[i] / 2 * 2 + 1;
            else {
                if (i == getDim() - 1)
                    result *= rSpaceSize[i] / 2 + 1;
                else
                    result *= rSpaceSize[i] / 2 * 2 + 1;
            }
                
        }
        return result;
    }

    template<class ScalarType, size_t Dim>
    Utils::Array<ssize_t, Dim> FFT<ScalarType, Dim>::linearIndexToDim(size_t index) const {
        Utils::Array<ssize_t, Dim> result(getDim());
        for (size_t i = 0; i < getDim(); ++i) {
            const size_t componentsSizeFrom_i = componentsSizeFrom(i + 1);
            ssize_t dim_i = index / componentsSizeFrom_i;
            index -= componentsSizeFrom_i * dim_i;
            if constexpr (isComplex)
                result[i] = dim_i - rSpaceSize[i] / 2;
            else {
                if (i == getDim() - 1)
                    result[i] = dim_i;
                else
                    result[i] = dim_i - rSpaceSize[i] / 2;
            }
        }
        return result;
    }

    template<class ScalarType, size_t Dim>
    bool FFT<ScalarType, Dim>::checkSize(const Utils::Array<size_t, Dim>& rSpaceSize) {
        for (size_t elem : rSpaceSize)
            if (elem > static_cast<size_t>(std::numeric_limits<int>::max()))
                return false;
        return true;
    }
}
