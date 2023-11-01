include(CTest)
include(../cmake/Modules/PhysicaAddTest.cmake)

message(STATUS "Configuring tests.")

# Core
## AI
physica_add_test(AP Core/AI/AP.cpp)
## Interface
physica_add_test(VaspWarpper Core/Interface/VaspWarpper.cpp)
## IO
physica_add_test(Gnuplot Core/IO/Gnuplot.cpp)
physica_add_test(HDF5 Core/IO/HDF5.cpp)
physica_add_test(Poscar Core/IO/Poscar.cpp)
physica_add_test(PWscfIn Core/IO/PWscfIn.cpp)
## Math
### Algebra
#### Equation Solver
physica_add_test(ElementaryEquation Core/Math/Algebra/EquationSolver/ElementaryEquation.cpp)
#### Linear Algebra
physica_add_test(EigenSolver Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.cpp)
physica_add_test(JacobiDavidson Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.cpp)
physica_add_test(SymmEigenSolver Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.cpp)
physica_add_test(Bidiagonalization Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Bidiagonalization.cpp)
physica_add_test(Cholesky Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/CholeskyTest.cpp)
physica_add_test(Hessenburg Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/HessenburgTest.cpp)
physica_add_test(LUDecomposition Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/LUDecomposition.cpp)
physica_add_test(Schur Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Schur.cpp)
physica_add_test(SVD Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/SVD.cpp)
physica_add_test(Tridiagonalization Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.cpp)
physica_add_test(DenseMatrix Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cpp)
physica_add_test(DenseMatrixExpression Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixExpression.cpp)
physica_add_test(InverseMatrix Core/Math/Algebra/LinearAlgebra/Matrix/InverseMatrixTest.cpp)
physica_add_test(MatrixProduct Core/Math/Algebra/LinearAlgebra/Matrix/MatrixProduct.cpp)
physica_add_test(PermutationMatrix Core/Math/Algebra/LinearAlgebra/Matrix/PermutationMatrix.cpp)
physica_add_test(SparseMatrix Core/Math/Algebra/LinearAlgebra/Matrix/SparseMatrix.cpp)
physica_add_test(Givens Core/Math/Algebra/LinearAlgebra/Vector/Givens.cpp)
physica_add_test(HouseHolder Core/Math/Algebra/LinearAlgebra/Vector/HouseHolder.cpp)
physica_add_test(Vector Core/Math/Algebra/LinearAlgebra/Vector/Vector.cpp)
physica_add_test(LinearEquations Core/Math/Algebra/LinearAlgebra/LinearEquations/LinearEquations.cpp)
physica_add_test(LUSolver Core/Math/Algebra/LinearAlgebra/LinearEquations/LUSolver.cpp)
physica_add_test(Polynomial Core/Math/Algebra/Polynomial.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(DenseMatrix_cuda Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cu)
    physica_add_test(Vector_cuda Core/Math/Algebra/LinearAlgebra/Vector/Vector.cu)
endif()
### Calculus
physica_add_test(ODESolver Core/Math/Calculus/ODE/ODESolver.cpp)
physica_add_test(SRK2 Core/Math/Calculus/ODE/SRK2.cpp)
physica_add_test(PoissonFEM Core/Math/Calculus/PDE/FEM/PoissonFEM.cpp)
physica_add_test(Chebyshev Core/Math/Calculus/ChebyshevTest.cpp)
physica_add_test(Differential Core/Math/Calculus/Differential.cpp)
physica_add_test(Integrate Core/Math/Calculus/IntegrateTest.cpp)
physica_add_test(Interpolation Core/Math/Calculus/Interpolation.cpp)
physica_add_test(SpetialFunctions Core/Math/Calculus/SpetialFunctionsTest.cpp)
### Geometry
physica_add_test(GeoBase Core/Math/Geometry/GeoBase.cpp)
physica_add_test(CubeCross Core/Math/Geometry/CubeCross.cpp)
### Optimization
physica_add_test(ConjugateGradient Core/Math/Optimization/ConjugateGradient.cpp)
physica_add_test(EqualityQuadraticProgramming Core/Math/Optimization/EqualityQuadraticProgramming.cpp)
physica_add_test(QuadraticProgramming Core/Math/Optimization/QuadraticProgramming.cpp)
physica_add_test(SteepestDescent Core/Math/Optimization/SteepestDescent.cpp)
### Statistics
physica_add_test(NumCharacter Core/Math/Statistics/NumCharacter.cpp)
### Transform
physica_add_test(FFT Core/Math/Transform/FFT.cpp)
## MultiPrecision
physica_add_test(Const Core/MultiPrecision/ConstTest.cpp)
physica_add_test(Differentiable Core/MultiPrecision/Differentiable.cpp)
physica_add_test(Integer Core/MultiPrecision/IntegerTest.cpp)
physica_add_test(Scalar Core/MultiPrecision/ScalarTest.cpp)
## Parallel
physica_add_test(ThreadPool Core/Parallel/ThreadPool.cpp)
## Physics
physica_add_test(GaussBase Core/Physics/GaussBase.cpp)
physica_add_test(IceGenerator Core/Physics/IceGenerator.cpp)
physica_add_test(RHFSolver Core/Physics/RHFSolver.cpp)
### MD
physica_add_test(PairModel Core/Physics/MD/ForceModel/PairModel.cpp)
physica_add_test(SilveraGoldman Core/Physics/MD/ForceModel/SilveraGoldman.cpp)
physica_add_test(TodaModel Core/Physics/MD/ForceModel/TodaModel.cpp)
physica_add_test(Langevin Core/Physics/MD/Thermostat/Langevin.cpp)
physica_add_test(MDCell Core/Physics/MD/MDCell.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(HardCore_cuda Core/Physics/MD/HardCore.cu)
endif()
### SolidState
physica_add_test(BrillouInterpolate Core/Physics/SolidState/BrillouInterpolate.cpp)
physica_add_test(CrystalCell Core/Physics/SolidState/CrystalCell.cpp)
# Logger
physica_add_test(Logger Logger/LoggerTest.cpp)
# Utils
physica_add_test(DirStack Utils/DirStackTest.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
