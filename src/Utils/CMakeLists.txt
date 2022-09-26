##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_UTILS *.cpp)
if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_UTILS *.cu)
    list(APPEND SOURCES_UTILS ${CUDA_SOURCES_UTILS})
endif()
add_library(PhysicaUtils SHARED ${SOURCES_UTILS})

install(TARGETS PhysicaUtils DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Utils DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaUtils)
