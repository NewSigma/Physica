# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE DeviceProp.cu)
endif()
