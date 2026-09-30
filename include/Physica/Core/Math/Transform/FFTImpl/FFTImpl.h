/*
 * Copyright 2020-2025 Weibo He.
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

#include "ThreadGuardFFTW.h"
#include "../FFT.h"

namespace Physica {
    template<Scalar T>
    FFT<T, 1>::FFT()
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(nullptr)
            , rSpaceSize(0)
            , planFlag(PlanFlag::Measure) {}

    template<Scalar T>
    FFT<T, 1>::FFT(size_t rSpaceSize_)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , rSpaceSize(static_cast<int>(rSpaceSize_))
            , planFlag(PlanFlag::Measure) {
        assert(rSpaceSize_ <= static_cast<size_t>(std::numeric_limits<int>::max()));
        buffer = static_cast<ComplexType*>(Internal::fftw_malloc(getKSpaceSize() * sizeof(ComplexType)));
    }

    template<Scalar T>
    FFT<T, 1>::FFT(size_t rSpaceSize_, PlanFlag planFlag_)
            : FFT(rSpaceSize_) {
        planFlag = planFlag_;
        initializePlan();
    }

    template<Scalar T>
    FFT<T, 1>::FFT(const VectorND<ScalarType>& data_, PlanFlag planFlag)
            : FFT(data_.getLength(), planFlag) {
        transform(data_);
    }

    template<Scalar T>
    FFT<T, 1>::FFT(const FFT& fft)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(static_cast<ComplexType*>(Internal::fftw_malloc(fft.rSpaceSize * sizeof(ComplexType))))
            , rSpaceSize(fft.rSpaceSize)
            , planFlag(fft.planFlag) {
        initializePlan();
    }

    template<Scalar T>
    FFT<T, 1>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(fft.rSpaceSize)
            , planFlag(fft.planFlag) {
        fft.forward_plan = FFTPlan(nullptr);
        fft.backward_plan = FFTPlan(nullptr);
        fft.buffer = nullptr;
    }

    template<Scalar T>
    FFT<T, 1>::~FFT() noexcept {
        std::unique_lock<std::mutex> locker(Internal::ThreadGuardFFTW::getInstance().globalMutex);
        Internal::fftw_destroy_plan<isSinglePrec()>(forward_plan);
        Internal::fftw_destroy_plan<isSinglePrec()>(backward_plan);
        Internal::fftw_free(buffer);
    }

    template<Scalar T>
    void FFT<T, 1>::swap(FFT& __restrict fft) noexcept {
        assert(this != &fft && "[Error]: Self swap is likely a bug");
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        std::swap(rSpaceSize, fft.rSpaceSize);
        std::swap(planFlag, fft.planFlag);
    }

    template<Scalar T>
    FFT<T, 1>::RealType FFT<T, 1>::getRSpaceDelta(RealType kSpaceDelta) const noexcept {
        return RealType(2 * M_PI) / (kSpaceDelta * getRSpaceSize());
    }

    template<Scalar T>
    FFT<T, 1>::RealType FFT<T, 1>::getKSpaceDelta(RealType rSpaceDelta) const noexcept {
        return RealType(2 * M_PI) / (rSpaceDelta * getRSpaceSize());
    }

    template<Scalar T>
    __host__ __device__ consteval bool FFT<T, 1>::isComplex() noexcept {
        return Traits<This>::isComplex;
    }

    template<Scalar T>
    __host__ __device__ consteval bool FFT<T, 1>::isSinglePrec() noexcept {
        return Traits<This>::isSinglePrec;
    }

    template<Scalar T>
    FFT<T, 1> FFT<T, 1>::makeEmptyFFT(size_t rSpaceSize) {
        return FFT<T, 1>(rSpaceSize);
    }

    template<Scalar T>
    template<std::integral IndexType>
    __host__ __device__ constexpr IndexType FFT<T, 1>::rSizeToKSize(IndexType rSize) noexcept {
        if constexpr (isComplex())
            return rSize;
        else
            return rSize / 2 + 1;
    }

    template<Scalar T>
    void FFT<T, 1>::transform(const This& planProvider, This& bufferProvider) {
        const auto forward_plan = planProvider.forward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(forward_plan != FFTPlan(nullptr) && "[Error]: Bad plan provider or working on a empry fft");
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isComplex())
            Internal::fftw_execute_dft<isSinglePrec()>(forward_plan, buffer);
        else
            Internal::fftw_execute_dft_r2c<isSinglePrec()>(forward_plan, buffer);
    }

    template<Scalar T>
    void FFT<T, 1>::rawInvTransform(const This& planProvider, This& bufferProvider) {
        const auto backward_plan = planProvider.backward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(backward_plan != FFTPlan(nullptr));
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isComplex())
            Internal::fftw_execute_dft<isSinglePrec()>(backward_plan, buffer);
        else
            Internal::fftw_execute_dft_c2r<isSinglePrec()>(backward_plan, buffer);
    }

    template<Scalar T>
    void FFT<T, 1>::invTransform(const This& planProvider, This& bufferProvider) {
        rawInvTransform(planProvider, bufferProvider);
        const ScalarType factor = RealType(1.0 / planProvider.getRSpaceSize());
        bufferProvider.getRSpace() *= factor;
    }

    template<Scalar T>
    void FFT<T, 1>::initializePlan() noexcept {
        assert(forward_plan == FFTPlan(nullptr));
        assert(backward_plan == FFTPlan(nullptr));
        std::unique_lock<std::mutex> locker(Internal::ThreadGuardFFTW::getInstance().globalMutex);
        const auto dims = Array<int, 1>::generate([this](size_t) { return rSpaceSize; });
        if constexpr (isComplex()) {
            forward_plan = Internal::fftw_plan_dft<isSinglePrec()>(1, dims.data(), buffer, true, planFlag);
            backward_plan = Internal::fftw_plan_dft<isSinglePrec()>(1, dims.data(), buffer, false, planFlag);
        }
        else {
            forward_plan = Internal::fftw_plan_dft_r2c<isSinglePrec()>(1, dims.data(), buffer, planFlag);
            backward_plan = Internal::fftw_plan_dft_c2r<isSinglePrec()>(1, dims.data(), buffer, planFlag);
        }
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::FFT()
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(nullptr)
            , rSpaceSize(Dim, 0)
            , kSpaceSize(Dim, 0)
            , planFlag(PlanFlag::Measure) {}

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::FFT(const Array<size_t, Dim>& rSpaceSize_)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , rSpaceSize(rSpaceSize_.getLength())
            , planFlag(PlanFlag::Measure) {
        assert(checkSize(rSpaceSize_));
        for (size_t i = 0; i < rSpaceSize_.getLength(); ++i)
            rSpaceSize[i] = static_cast<int>(rSpaceSize_[i]);
        kSpaceSize = rSizeToKSize(rSpaceSize);

        buffer = static_cast<ComplexType*>(Internal::fftw_malloc(sumKSpaceSize(0) * sizeof(ComplexType)));
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::FFT(const Array<size_t, Dim>& rSpaceSize_, PlanFlag planFlag_)
            : FFT(rSpaceSize_) {
        planFlag = planFlag_;
        initializePlan();
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::FFT(const FFT& fft)
            : forward_plan(nullptr)
            , backward_plan(nullptr)
            , buffer(static_cast<ComplexType*>(Internal::fftw_malloc(fft.sumKSpaceSize(0) * sizeof(ComplexType))))
            , rSpaceSize(fft.rSpaceSize)
            , kSpaceSize(fft.kSpaceSize)
            , planFlag(fft.planFlag) {
        initializePlan();
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::FFT(FFT&& fft) noexcept
            : forward_plan(fft.forward_plan)
            , backward_plan(fft.backward_plan)
            , buffer(fft.buffer)
            , rSpaceSize(std::move(fft.rSpaceSize))
            , kSpaceSize(std::move(fft.kSpaceSize))
            , planFlag(fft.planFlag) {
        fft.forward_plan = FFTPlan(nullptr);
        fft.backward_plan = FFTPlan(nullptr);
        fft.buffer = nullptr;
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::~FFT() noexcept {
        std::unique_lock<std::mutex> locker(Internal::ThreadGuardFFTW::getInstance().globalMutex);
        Internal::fftw_destroy_plan<isSinglePrec()>(forward_plan);
        Internal::fftw_destroy_plan<isSinglePrec()>(backward_plan);
        Internal::fftw_free(buffer);
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::swap(FFT& __restrict fft) noexcept {
        assert(this != &fft && "[Error]: Self swap is likely a bug");
        std::swap(forward_plan, fft.forward_plan);
        std::swap(backward_plan, fft.backward_plan);
        std::swap(buffer, fft.buffer);
        rSpaceSize.swap(fft.rSpaceSize);
        kSpaceSize.swap(fft.kSpaceSize);
        std::swap(planFlag, fft.planFlag);
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::RealType FFT<T, Dim>::getRSpaceDelta(
            RealType kSpaceDelta, unsigned int dim) const noexcept{
        return RealType(2 * M_PI) / (kSpaceDelta * getRSpaceSize()[dim]);
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim>::RealType FFT<T, Dim>::getKSpaceDelta(
            RealType rSpaceDelta, unsigned int dim) const noexcept {
        return RealType(2 * M_PI) / (rSpaceDelta * getRSpaceSize()[dim]);
    }

    template<Scalar T, size_t Dim>
    __host__ __device__ consteval bool FFT<T, Dim>::isComplex() noexcept {
        return Traits<This>::isComplex;
    }

    template<Scalar T, size_t Dim>
    __host__ __device__ consteval bool FFT<T, Dim>::isSinglePrec() noexcept {
        return Traits<This>::isSinglePrec;
    }

    template<Scalar T, size_t Dim>
    FFT<T, Dim> FFT<T, Dim>::makeEmptyFFT(const Array<size_t, Dim>& rSpaceSize) {
        return FFT(rSpaceSize);
    }

    template<Scalar T, size_t Dim>
    template<std::integral IndexType>
    Array<IndexType, Dim> FFT<T, Dim>::rSizeToKSize(const Array<IndexType, Dim>& rSize) {
        Array<IndexType, Dim> result(rSize.getLength());
        size_t i = 0;
        for (; i < rSize.getLength() - 1; ++i)
            result[i] = rSize[i];
        result[i] = FFT<T>::rSizeToKSize(rSize[i]);
        return result;
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::transform(const This& planProvider, This& bufferProvider) {
        const auto forward_plan = planProvider.forward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(forward_plan != FFTPlan(nullptr) && "[Error]: Bad plan provider or working on a empry fft");
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isComplex())
            Internal::fftw_execute_dft<isSinglePrec()>(forward_plan, buffer);
        else
            Internal::fftw_execute_dft_r2c<isSinglePrec()>(forward_plan, buffer);
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::rawInvTransform(const This& planProvider, This& bufferProvider) {
        const auto backward_plan = planProvider.backward_plan;
        const auto buffer = bufferProvider.buffer;
        assert(backward_plan != FFTPlan(nullptr));
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        if constexpr (isComplex())
            Internal::fftw_execute_dft<isSinglePrec()>(backward_plan, buffer);
        else
            Internal::fftw_execute_dft_c2r<isSinglePrec()>(backward_plan, buffer);
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::invTransform(const This& planProvider, This& bufferProvider) {
        rawInvTransform(planProvider, bufferProvider);
        const ScalarType factor = RealType(1.0 / planProvider.sumRSpaceSize(0));
        bufferProvider.getRSpace() *= factor;
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::initializePlan() noexcept {
        assert(forward_plan == FFTPlan(nullptr));
        assert(backward_plan == FFTPlan(nullptr));
        std::unique_lock<std::mutex> locker(Internal::ThreadGuardFFTW::getInstance().globalMutex);
        forward_plan = makeForwardPlan();
        backward_plan = makeBackwardPlan();
    }

    template<Scalar T, size_t Dim>
    FFTPlan FFT<T, Dim>::makeForwardPlan() {
        if constexpr (Dim == 2) {
            const auto dims = Array<int, 2>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); });
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(2, dims.data(), buffer, true, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_r2c<isSinglePrec()>(2, dims.data(), buffer, PlanFlag::Estimate);
        }
        else if constexpr (Dim == 3) {
            const auto dims = Array<int, 3>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); });
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(3, dims.data(), buffer, true, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_r2c<isSinglePrec()>(3, dims.data(), buffer, PlanFlag::Estimate);
        }
        else {
            const int rank = static_cast<int>(getDim());
            const auto dims = Array<int, Dynamic>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); }, getDim());
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(rank, dims.data(), buffer, true, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_r2c<isSinglePrec()>(rank, dims.data(), buffer, PlanFlag::Estimate);
        }
    }

    template<Scalar T, size_t Dim>
    FFTPlan FFT<T, Dim>::makeBackwardPlan() {
        if constexpr (Dim == 2) {
            const auto dims = Array<int, 2>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); });
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(2, dims.data(), buffer, false, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_c2r<isSinglePrec()>(2, dims.data(), buffer, PlanFlag::Estimate);
        }
        else if constexpr (Dim == 3) {
            const auto dims = Array<int, 3>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); });
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(3, dims.data(), buffer, false, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_c2r<isSinglePrec()>(3, dims.data(), buffer, PlanFlag::Estimate);
        }
        else {
            const int rank = static_cast<int>(getDim());
            const auto dims = Array<int, Dynamic>::generate([this](size_t i) { return static_cast<int>(rSpaceSize[i]); }, getDim());
            if constexpr (isComplex())
                return Internal::fftw_plan_dft<isSinglePrec()>(rank, dims.data(), buffer, false, PlanFlag::Estimate);
            else
                return Internal::fftw_plan_dft_c2r<isSinglePrec()>(rank, dims.data(), buffer, PlanFlag::Estimate);
        }
    }

    template<Scalar T, size_t Dim>
    size_t FFT<T, Dim>::sumRSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDim(); ++i)
            result *= getRSpaceSize()[i];
        return result;
    }

    template<Scalar T, size_t Dim>
    size_t FFT<T, Dim>::sumKSpaceSize(size_t from_dim) const {
        size_t result = 1;
        for (size_t i = from_dim; i < getDim(); ++i)
            result *= getKSpaceSize()[i];
        return result;
    }

    template<Scalar T, size_t Dim>
    void FFT<T, Dim>::normalizeIndexes(Array<ssize_t, Dim>& indexes) const {
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

    template<Scalar T, size_t Dim>
    size_t FFT<T, Dim>::componentsSizeFrom(size_t dim) const {
        size_t result = 1;
        for (size_t i = dim; i < getDim(); ++i) {
            if constexpr (isComplex())
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

    template<Scalar T, size_t Dim>
    Array<ssize_t, Dim> FFT<T, Dim>::linearIndexToDim(size_t index) const {
        Array<ssize_t, Dim> result(getDim());
        for (size_t i = 0; i < getDim(); ++i) {
            const size_t componentsSizeFrom_i = componentsSizeFrom(i + 1);
            ssize_t dim_i = index / componentsSizeFrom_i;
            index -= componentsSizeFrom_i * dim_i;
            if constexpr (isComplex())
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

    template<Scalar T, size_t Dim>
    bool FFT<T, Dim>::checkSize(const Array<size_t, Dim>& rSpaceSize) {
        return std::all_of(rSpaceSize.begin(), rSpaceSize.end(), [](size_t elem) {
            return elem <= static_cast<size_t>(std::numeric_limits<int>::max());
        });
    }
}
