/*
 * Copyright 2023 WeiBo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Random/RandomPool.h"

using namespace Physica::Core;
using ScalarType = ComplexScalar<Scalar<Double>>;
using MatrixType = DenseHermiteMatrix<ScalarType, Dynamic, Dynamic>;
using RandomPoolType = RandomPool<std::mt19937>;

int main() {
    MatrixType mat(3);
    mat.random_uniform(RandomPoolType::getInstance().getGen());
    std::cout << mat.asMatrix() << std::endl;
    return 0;
}
