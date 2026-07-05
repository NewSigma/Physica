if (${PHYSICA_HDF5})
    add_subdirectory(HDF5)
endif()

target_sources(PhysicaCore PRIVATE
    Gnuplot.cpp
    Mnist.cpp
)
