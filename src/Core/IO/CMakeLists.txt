if (${PHYSICA_HDF5})
    add_subdirectory(HDF5)
endif()
add_subdirectory(QE)
add_subdirectory(VASP)

target_sources(PhysicaCore PRIVATE
    Gnuplot.cpp
    Mnist.cpp
)
