# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025 Weibo He
# SPDX-License-Identifier: MPL-2.0
physica_add_example(DenseVector DenseVector.cpp)
physica_add_example(JacobiDavidson JacobiDavidson.cpp)

if(${PHYSICA_CUDA})
    physica_add_example(DenseVector_cuda DenseVector.cu)
endif()
