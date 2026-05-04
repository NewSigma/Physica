target_sources(PhysicaCore PRIVATE IntegerImpl/Math.cpp)

add_subdirectory(RealImpl)

target_sources(PhysicaCore PRIVATE ArrayArithmetic.cpp Integer.cpp Rational.cpp Real.cpp)
