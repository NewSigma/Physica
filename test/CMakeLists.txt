include(CTest)
include(../cmake/Modules/PhysicaAddTest.cmake)

message(STATUS "Configuring tests.")
include_directories(./)
# Math
## Algebra
### Equation Solver
physica_add_test(ElementaryEquation Core/Math/Algebra/EquationSolver/ElementaryEquation.cpp)
### Linear Algebra
physica_add_test(Cholesky Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/CholeskyTest.cpp)
physica_add_test(Hessenburg Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/HessenburgTest.cpp)
physica_add_test(RealSchur Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/RealSchurTest.cpp)
physica_add_test(DenseMatrixExpression Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixExpression.cpp)
physica_add_test(InverseMatrix Core/Math/Algebra/LinearAlgebra/Matrix/InverseMatrixTest.cpp)
physica_add_test(MatrixProduct Core/Math/Algebra/LinearAlgebra/Matrix/MatrixProduct.cpp)
physica_add_test(EigenSolver Core/Math/Algebra/LinearAlgebra/EigenSolver.cpp)
physica_add_test(Givens Core/Math/Algebra/LinearAlgebra/GivensTest.cpp)
physica_add_test(HouseHolder Core/Math/Algebra/LinearAlgebra/HouseHolderTest.cpp)
physica_add_test(LinearEquations Core/Math/Algebra/LinearAlgebra/LinearEquations.cpp)
physica_add_test(LUDecomposition Core/Math/Algebra/LinearAlgebra/LUDecompositionTest.cpp)
physica_add_test(Vector Core/Math/Algebra/LinearAlgebra/VectorTest.cpp)
physica_add_test(Polynomial Core/Math/Algebra/Polynomial.cpp)
## Calculus
physica_add_test(Chebyshev Core/Math/Calculus/ChebyshevTest.cpp)
physica_add_test(Differential Core/Math/Calculus/Differential.cpp)
physica_add_test(Function Core/Math/Calculus/FunctionTest.cpp)
physica_add_test(Integrate Core/Math/Calculus/IntegrateTest.cpp)
physica_add_test(ODE Core/Math/Calculus/ODETest.cpp)
physica_add_test(SpetialFunctions Core/Math/Calculus/SpetialFunctionsTest.cpp)
## Optimization
physica_add_test(Adam Core/Math/Optimization/AdamTest.cpp)
physica_add_test(ConjugateGradient Core/Math/Optimization/ConjugateGradientTest.cpp)
physica_add_test(EqualityQuadraticProgramming Core/Math/Optimization/EqualityQuadraticProgramming.cpp)
physica_add_test(QuadraticProgramming Core/Math/Optimization/QuadraticProgramming.cpp)
# MultiPrecision
physica_add_test(Const Core/MultiPrecision/ConstTest.cpp)
physica_add_test(Integer Core/MultiPrecision/IntegerTest.cpp)
physica_add_test(Scalar Core/MultiPrecision/ScalarTest.cpp)
# Parallel
physica_add_test(ThreadPool Core/Parallel/ThreadPoolTest.cpp)
# Physics
physica_add_test(GaussBase Core/Physics/GaussBase.cpp)
physica_add_test(RHFSolver Core/Physics/RHFSolver.cpp)
# Logger
physica_add_test(Logger Logger/LoggerTest.cpp)
# Utils
physica_add_test(DirStack Utils/DirStackTest.cpp)
