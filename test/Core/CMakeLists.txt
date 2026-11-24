# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Scalar)
## ML
physica_add_test(AP ML/AP.cpp)
physica_add_test(BayesOpt ML/BayesOpt.cpp)
physica_add_test(GPModel ML/GPModel.cpp)

add_subdirectory(Parallel)
## Utils
physica_add_test(Array Utils/Array.cpp)
physica_add_test(ArrayND Utils/ArrayND.cpp)
physica_add_test(Handle Utils/Handle.cpp)
physica_add_test(Range Utils/Range.cpp)

add_subdirectory(Physics)
