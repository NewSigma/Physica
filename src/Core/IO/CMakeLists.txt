if (${PHYSICA_HDF5})
    add_subdirectory(HDF5)
endif()
add_subdirectory(QE)
add_subdirectory(VASP)

target_sources(PhysicaCore PRIVATE
    CsvFile.cpp
    Gnuplot.cpp
    Mnist.cpp
)
