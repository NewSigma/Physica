/*
 * Copyright 2022 WeiBo He.
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

#include <torch/torch.h>
#include "Physica/Utils/Template/CRTPBase.h"
#include "HyperParam/ParamSet.h"

namespace Physica::AI {
    namespace Internal {
        template<class T> class Traits;
    }

    template<class Derived>
    class Model : public Utils::CRTPBase<Derived>, public torch::nn::Module {
    public:
        using Base = torch::nn::Module;
        using DataSet = typename Internal::Traits<Derived>::DataSet;
        using HyperParams = typename Internal::Traits<Derived>::HyperParams;

        HyperParams active_params;
    public:
        /* Operations */
        void init() { this->getDerived().init(); }
        void train(const DataSet& data) { this->getDerived().train(data); }
        /* Getters */
        torch::Tensor forward(torch::Tensor x) { return this->getDerived().forward(std::move(x)); }
        double loss(const torch::Tensor& features, const torch::Tensor& labels) { return this->getDerived().loss(features, labels); }
    protected:
        Model() = default;
    };
}
