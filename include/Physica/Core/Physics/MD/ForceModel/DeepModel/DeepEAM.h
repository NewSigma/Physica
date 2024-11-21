/*
 * Copyright 2024 Weibo He.
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

#include "Descriptor/ChebyshevRadial.h"
#include "ConservedFieldNet.h"

namespace Physica::Core {
    template<class NetType, bool IsSmallCell>
    class DeepEAM {
        static_assert(std::is_base_of<ConservedFieldNet<NetType>, NetType>::value, "[Error]: DeepEAM must be a conserved field");
        using ScalarType = NetType::ScalarType;
        using ValueType = ScalarType::ValueType;
        using DescriptorBase = ChebyshevRadial<ScalarType, IsSmallCell>;
        using MDCellType = DescriptorBase::MDCellType;

        DescriptorBase base;
        NetType net;
    public:
        DeepEAM(DescriptorBase base_, NetType net_);
        DeepEAM(const DeepEAM&) = default;
        DeepEAM(DeepEAM&&) noexcept = default;
        ~DeepEAM() = default;
        /* Operators */
        DeepEAM& operator=(DeepEAM obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ValueType potentialV(const MDCellType& cell) const;

        template<class Executor>
        [[nodiscard]] VectorND<ValueType> force(const MDCellType& cell) const;
        template<Vector V, class Executor>
        void forceAsync([[maybe_unused]] const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        void swap(DeepEAM& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const NetType& getNet() const noexcept { return net; }
    };

    template<class NetType, bool IsSmallCell>
    DeepEAM<NetType, IsSmallCell>::DeepEAM(DescriptorBase base_, NetType net_) : base(std::move(base_)), net(std::move(net_)) {}

    template<class NetType, bool IsSmallCell>
    DeepEAM<NetType, IsSmallCell>::ValueType
    DeepEAM<NetType, IsSmallCell>::potentialV(const MDCellType& cell) const {
        ValueType result = 0;
        const auto descriptors = base.project(cell);
        for (size_t i = 0; i < cell.getNumParticle(); ++i)
            result += net.forward(descriptors[i].flatten()).getValue();
        return result;
    }

    template<class NetType, bool IsSmallCell>
    void DeepEAM<NetType, IsSmallCell>::swap(DeepEAM& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        base.swap(obj.base);
        net.swap(obj.net);
    }
}
