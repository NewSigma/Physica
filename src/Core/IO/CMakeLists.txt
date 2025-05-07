if (${PHYSICA_HDF5})
    file(GLOB_RECURSE SOURCES HDF5/*.cpp)
    target_sources(PhysicaCore PRIVATE ${SOURCES})
endif()

file(GLOB_RECURSE SOURCES QE/*.cpp VASP/*.cpp)
target_sources(PhysicaCore PRIVATE ${SOURCES} CsvFile.cpp Gnuplot.cpp Mnist.cpp)
