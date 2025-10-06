target_sources(PhysicaCore PRIVATE IntegerImpl/Math.cpp)

add_subdirectory(RealImpl)

target_sources(PhysicaCore PRIVATE Integer.cpp)
target_sources(PhysicaCore PRIVATE Rational.cpp)
