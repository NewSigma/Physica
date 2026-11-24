/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "LValueTensor.h"

namespace Physica {
    template<class Derived>
    class StridedTensor : public LValueTensor<Derived> {
        using Base = LValueTensor<Derived>;
        using This = StridedTensor<Derived>;
    public:
        using Base::NDim;
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    public:
        ~StridedTensor() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;

        [[nodiscard]] bool operator==(this const auto& self, const Tensor auto& other) noexcept;
        /* Operations */
        void assign(Tensor auto&& target) const noexcept;
        void zeros() noexcept;
        /* Getters */
        [[nodiscard]] auto getStrides() const noexcept;
        [[nodiscard]] size_t getStride(size_t dim) const noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, std::same_as<size_t> auto... dims) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile(size_t dim) noexcept;
    protected:
        StridedTensor() = default;
        StridedTensor(const This&) = default;
        StridedTensor(This&&) noexcept = default;
    private:
        /* Getters */
        template<class Self, class Other>
        [[nodiscard]] __host__ __device__ consteval static bool isMinorCompact() noexcept;
        [[nodiscard]] __host__ __device__ constexpr static int getCompactMinorDim(const auto& tensor) noexcept;
        [[nodiscard]] __host__ __device__ constexpr static size_t getCompactMinorSize(const auto& tensor) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static auto makeCompactBlock(auto&& tensor, const IndexType& index, int minorDim) noexcept;
        /* Friends */
        friend class device_obj<This>;
    };
}

#include "StridedTensorImpl/StridedTensorImpl.h"
#include "StridedTensorImpl/TensorFiber.h"
