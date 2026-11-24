# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
physica_add_example(DecisionTree DecisionTree.cpp)
physica_add_example(GPModel GPModel.cpp)
physica_add_example(Mnist Mnist.cpp)
if(${PHYSICA_CUDA})
    physica_add_example(Mnist_cuda Mnist.cu)
endif()
