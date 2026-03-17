##############################################Libs################################################
find_package(LLVM REQUIRED CONFIG)

include(AddLLVM)
#############################################Project##############################################
add_llvm_pass_plugin(PhysicaTransforms
  Plugin.cpp
  LifetimeMovePass.cpp

  LINK_COMPONENTS
  Analysis
  Core
  ProfileData
  Support
  TargetParser
)

target_include_directories(PhysicaTransforms SYSTEM PRIVATE ${LLVM_INCLUDE_DIRS})
