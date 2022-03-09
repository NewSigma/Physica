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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/LValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::AI {
    template<class VectorType>
    torch::Tensor toTensor(const Physica::Core::LValueVector<VectorType>& v, at::TensorOptions options = {}) {
        const int64_t length = static_cast<int64_t>(v.getLength());
        torch::Tensor result = torch::empty({length}, std::move(options));
        for (int64_t i = 0; i < length; ++i)
            result[i] = v[i].getTrivial();
        return torch::Tensor(std::move(result));
    }

    template<class MatrixType>
    torch::Tensor toTensor(const Physica::Core::LValueMatrix<MatrixType>& m, at::TensorOptions options = {}) {
        const int64_t row = static_cast<int64_t>(m.getRow());
        const int64_t column = static_cast<int64_t>(m.getColumn());
        torch::Tensor result = torch::empty({row, column}, std::move(options));
        for (int64_t i = 0; i < row; ++i)
            for (int64_t j = 0; j < column; ++j)
                result[i][j] = m(i, j).getTrivial();
        return torch::Tensor(std::move(result));
    }
}
