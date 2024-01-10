##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_PYTHON *.cpp)
if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_PYTHON *.cu)
    list(APPEND SOURCES_PYTHON ${CUDA_SOURCES_PYTHON})
endif()
pybind11_add_module(PhysicaPython SHARED ${SOURCES_PYTHON})

install(TARGETS PhysicaPython DESTINATION lib)
#############################################Options##############################################
target_link_libraries(PhysicaPython PRIVATE PhysicaCore)
if(${PHYSICA_GUI})
    target_link_libraries(PhysicaPython PRIVATE PhysicaGui)
endif()
