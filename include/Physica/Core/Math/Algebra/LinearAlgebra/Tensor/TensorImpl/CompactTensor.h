/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "StridedTensor.h"

namespace Physica {
    template<class Derived>
    class CompactTensor : public StridedTensor<Derived> {
        using Base = StridedTensor<Derived>;
        using This = CompactTensor<Derived>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    public:
        ~CompactTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;

        [[nodiscard]] bool operator==(this const auto& self, const Tensor auto& other) noexcept;
        /* Operations */
        void assign(Tensor auto&& target) const noexcept;
        /* Getters */
        [[nodiscard]] auto data(this auto&& self) noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] constexpr auto getStrides() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept { return true; }
    protected:
        CompactTensor() = default;
        CompactTensor(const This&) = default;
        CompactTensor(This&&) noexcept = default;
    private:
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static IndexType getStridesImpl(const IndexType& shape) noexcept;
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<CompactTensor<Derived>> : public Traits<Derived> {};
}

#include "CompactTensorImpl/CompactTensorImpl.h"
#include "CompactTensorImpl/Flatten.h"
