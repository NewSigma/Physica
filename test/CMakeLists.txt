include(CTest)
include(../cmake/Modules/PhysicaAddTest.cmake)

message(STATUS "Configuring tests.")
include_directories(./)
# IO
physica_add_test(Poscar Core/IO/Poscar.cpp)
physica_add_test(Gnuplot Core/IO/Gnuplot.cpp)
# Math
## Algebra
### Equation Solver
physica_add_test(ElementaryEquation Core/Math/Algebra/EquationSolver/ElementaryEquation.cpp)
### Linear Algebra
physica_add_test(Bidiagonalization Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Bidiagonalization.cpp)
physica_add_test(Cholesky Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/CholeskyTest.cpp)
physica_add_test(Hessenburg Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/HessenburgTest.cpp)
physica_add_test(LUDecomposition Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/LUDecompositionTest.cpp)
physica_add_test(Schur Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/SchurTest.cpp)
physica_add_test(Tridiagonalization Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.cpp)
physica_add_test(DenseMatrixExpression Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixExpression.cpp)
physica_add_test(InverseMatrix Core/Math/Algebra/LinearAlgebra/Matrix/InverseMatrixTest.cpp)
physica_add_test(MatrixProduct Core/Math/Algebra/LinearAlgebra/Matrix/MatrixProduct.cpp)
physica_add_test(Givens Core/Math/Algebra/LinearAlgebra/Vector/GivensTest.cpp)
physica_add_test(HouseHolder Core/Math/Algebra/LinearAlgebra/Vector/HouseHolderTest.cpp)
physica_add_test(Vector Core/Math/Algebra/LinearAlgebra/Vector/VectorTest.cpp)
physica_add_test(EigenSolver Core/Math/Algebra/LinearAlgebra/EigenSolver.cpp)
physica_add_test(LinearEquations Core/Math/Algebra/LinearAlgebra/LinearEquations.cpp)
physica_add_test(SymmEigenSolver Core/Math/Algebra/LinearAlgebra/SymmEigenSolver.cpp)
physica_add_test(Polynomial Core/Math/Algebra/Polynomial.cpp)
## Calculus
physica_add_test(ODE Core/Math/Calculus/ODE/ODETest.cpp)
physica_add_test(PoissonFEM Core/Math/Calculus/PDE/FEM/PoissonFEM.cpp)
physica_add_test(SRK2 Core/Math/Calculus/ODE/SRK2.cpp)
physica_add_test(Chebyshev Core/Math/Calculus/ChebyshevTest.cpp)
physica_add_test(Differential Core/Math/Calculus/Differential.cpp)
physica_add_test(Function Core/Math/Calculus/FunctionTest.cpp)
physica_add_test(Integrate Core/Math/Calculus/IntegrateTest.cpp)
physica_add_test(Interpolation Core/Math/Calculus/Interpolation.cpp)
physica_add_test(SpetialFunctions Core/Math/Calculus/SpetialFunctionsTest.cpp)
## Optimization
physica_add_test(ConjugateGradient Core/Math/Optimization/ConjugateGradientTest.cpp)
physica_add_test(EqualityQuadraticProgramming Core/Math/Optimization/EqualityQuadraticProgramming.cpp)
physica_add_test(QuadraticProgramming Core/Math/Optimization/QuadraticProgramming.cpp)
## Transform
physica_add_test(FFT Core/Math/Transform/FFT.cpp)
# MultiPrecision
physica_add_test(Const Core/MultiPrecision/ConstTest.cpp)
physica_add_test(Integer Core/MultiPrecision/IntegerTest.cpp)
physica_add_test(Scalar Core/MultiPrecision/ScalarTest.cpp)
# Parallel
physica_add_test(ThreadPool Core/Parallel/ThreadPoolTest.cpp)
# Physics
physica_add_test(Ewald Core/Physics/Ewald.cpp)
physica_add_test(GaussBase Core/Physics/GaussBase.cpp)
physica_add_test(RHFSolver Core/Physics/RHFSolver.cpp)
physica_add_test(PWBaseWave Core/Physics/PWBaseWave.cpp)
# Logger
physica_add_test(Logger Logger/LoggerTest.cpp)
# Utils
physica_add_test(DirStack Utils/DirStackTest.cpp)

if(NOT ${DEBUG_MODE})
    include(SlowTests.cmake)
endif()
