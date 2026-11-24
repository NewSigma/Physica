# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
##############################################Libs################################################
if(NOT PHYSICA_LLVM)
    message(FATAL_ERROR "Transforms requires LLVM package")
endif()

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
