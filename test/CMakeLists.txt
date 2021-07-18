include(CTest)
include(../cmake/Modules/PhysicaAddTest.cmake)

message(STATUS "Configuring tests.")
# Math
physica_add_test(Cholesky Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/CholeskyTest.cpp)
physica_add_test(HouseHolder Core/Math/Algebra/LinearAlgebra/HouseHolderTest.cpp)
physica_add_test(InverseMatrix Core/Math/Algebra/LinearAlgebra/InverseMatrixTest.cpp)
physica_add_test(LUDecomposition Core/Math/Algebra/LinearAlgebra/LUDecompositionTest.cpp)
physica_add_test(Vector Core/Math/Algebra/LinearAlgebra/VectorTest.cpp)
physica_add_test(Chebyshev Core/Math/Calculus/ChebyshevTest.cpp)
physica_add_test(Function Core/Math/Calculus/FunctionTest.cpp)
physica_add_test(Integrate Core/Math/Calculus/IntegrateTest.cpp)
physica_add_test(ODE Core/Math/Calculus/ODETest.cpp)
physica_add_test(SpetialFunctions Core/Math/Calculus/SpetialFunctionsTest.cpp)
physica_add_test(ConjugateGradient Core/Math/Optimization/ConjugateGradientTest.cpp)
# MultiPrecision
physica_add_test(ScalarAddSubExpression Core/MultiPrecision/ScalarImpl/ScalarAddSubExpressionTest.cpp)
physica_add_test(Const Core/MultiPrecision/ConstTest.cpp)
physica_add_test(Integer Core/MultiPrecision/IntegerTest.cpp)
physica_add_test(Scalar Core/MultiPrecision/ScalarTest.cpp)
# Parallel
physica_add_test(ThreadPool Core/Parallel/ThreadPoolTest.cpp)
# Logger
physica_add_test(Logger Logger/LoggerTest.cpp)
# Utils
physica_add_test(DirStack Utils/DirStackTest.cpp)
