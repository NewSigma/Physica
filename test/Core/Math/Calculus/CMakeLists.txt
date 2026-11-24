# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
# Integrate
physica_add_test(Vegas Integrate/Vegas.cpp)
physica_add_test(Integrate Integrate/Integrate.cpp)
# ODE
physica_add_test(ODESolver ODE/ODESolver.cpp)
physica_add_test(SRK2 ODE/SRK2.cpp)
# PDE
physica_add_test(PoissonFEM PDE/FEM/PoissonFEM.cpp)

physica_add_test(Chebyshev ChebyshevTest.cpp)
physica_add_test(Differential Differential.cpp)
physica_add_test(ElementaryEquation ElementaryEquation.cpp)
physica_add_test(Interpolation Interpolation.cpp)
physica_add_test(SpecialFunctions SpecialFunctionsTest.cpp)
