add_subdirectory(CUDA)
add_subdirectory(Unix)

target_sources(PhysicaCore PRIVATE
    BitArray.cpp
    Cycler.cpp
    NoImpl.cpp
    RingBuffer.cpp
    SSHWrapper.cpp
)
