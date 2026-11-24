# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2025 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(ML)
add_subdirectory(Math)
add_subdirectory(Physics)

physica_add_example(Diff Diff.cpp)
