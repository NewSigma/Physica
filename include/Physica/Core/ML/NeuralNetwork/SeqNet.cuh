/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "SeqNet.h"

namespace Physica {
    template<class Derived>
    class device_obj<SeqNet<Derived>> : public device_obj<LayerMixin<Derived>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: device_obj<> is unnecessary");
        using host_obj = SeqNet<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerMixin<Derived>>;
    public:
        using typename Base::ScalarType;
    protected:
        using Tv = ScalarType::ValueType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<RNG R, ExecutePolicy P>
        void train_step(int batchSize, const auto& dataset);
        template<RNG R, ExecutePolicy P>
        void train_step_for(int64_t numStep, int batchSize, const auto& dataset);

        [[nodiscard]] auto loss(const auto& dataset, size_t index) const { return Base::getDerived().loss(dataset, index); }
        [[nodiscard]] ScalarType loss(const auto& dataset) const;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    template<RNG R, ExecutePolicy P>
    void device_obj<SeqNet<Derived>>::train_step(int batchSize, const auto& dataset) {
        static_assert(Base::isTraining(), "[Error]: train_step must be called under training mode");
        const auto indices = R::random_int(batchSize, 0, dataset.getSize() - 1);
        const Tv mean_grad = reciprocal(Tv(batchSize));
        for (auto index : indices)
            loss(dataset, index).reverse(mean_grad);
        Base::step();
        Base::zero_grad();
    }

    template<class Derived>
    template<RNG R, ExecutePolicy P>
    void device_obj<SeqNet<Derived>>::train_step_for(int64_t numStep, int batchSize, const auto& dataset) {
        for (int64_t _ = 0; _ < numStep; ++_)
            train_step<R, P>(batchSize, dataset);
    }

    template<class Derived>
    auto device_obj<SeqNet<Derived>>::loss(const auto& dataset) const -> ScalarType {
        static_assert(Base::isInfering(), "[Error]: It is suggested using eval mode to reduce memory use");
        const size_t size = dataset.getSize();
        ScalarType result = 0;
        for (size_t i = 0; i < size; ++i)
            result.toNextMean(i, loss(dataset, i));
        return result;
    }
}
