physica_add_test(ThreadPool ThreadPool.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(ThreadBlock ThreadBlock.cu)
endif()
