add_subdirectory(Algebra)
add_subdirectory(Calculus)
# Discrete
physica_add_test(Combination Discrete/Combination.cpp)
# Geometry
physica_add_test(GeoBase Geometry/GeoBase.cpp)
physica_add_test(CubeCross Geometry/CubeCross.cpp)
# NumberTheory
physica_add_test(NumberTheory NumberTheory/NumberTheory.cpp)
# Optimization
physica_add_test(ConjugateGradient Optimization/ConjugateGradient.cpp)
physica_add_test(EqualityQuadraticProgramming Optimization/EqualityQuadraticProgramming.cpp)
physica_add_test(QuadraticProgramming Optimization/QuadraticProgramming.cpp)
physica_add_test(SteepestDescent Optimization/SteepestDescent.cpp)
# Statistics
physica_add_test(NumCharacter Statistics/NumCharacter.cpp)
physica_add_test(LinearFit Statistics/LinearFit.cpp)
# Transform
physica_add_test(FFT Transform/FFT.cpp)
