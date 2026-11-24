# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
target_sources(Benchmark PRIVATE DQMC.cpp)
target_sources(Benchmark PRIVATE ElasticDQMC.cpp)
target_sources(Benchmark PRIVATE FreqDQMC.cpp)
target_sources(Benchmark PRIVATE HubbardMatrix1D.cpp)
target_sources(Benchmark PRIVATE TPQ.cpp)
