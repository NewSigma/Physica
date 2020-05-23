include(Check/CheckIntegerCast.cmake)
include(Check/CheckEndian.cmake)
include(Check/CheckRightShift.cmake)
include(Check/CheckWordSize.cmake)
set(UseASM 1)
message(STATUS)

configure_file(PhysicaConfig.h.in ../PhysicaConfig.h)