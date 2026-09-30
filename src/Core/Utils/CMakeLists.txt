add_subdirectory(CUDA)
add_subdirectory(Unix)

if (${PHYSICA_MIMALLOC})
    target_sources(PhysicaCore PRIVATE Allocator/Mimalloc.cpp)
endif()

target_sources(PhysicaCore PRIVATE
    BitArray.cpp
    Cycler.cpp
    NoImpl.cpp
    RingBuffer.cpp
    Suspend.cpp
)
