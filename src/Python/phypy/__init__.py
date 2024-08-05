"""
Copyright 2024 Weibo He.

This file is part of Physica.

Physica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Physica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
"""
import sys
import os

sys.path.append(os.path.abspath('../../build/src/Python'))
"""
LLVM-17 uses dlsym() to locate symbols, and it is necessary import the symbols to python using os.RTLD_GLOBAL

Reference:
[1] https://github.com/apache/arrow/issues/39695
"""
def importPhysica():
    dlflags = sys.getdlopenflags()
    sys.setdlopenflags(dlflags + os.RTLD_GLOBAL)
    try:
        import PhysicaPython as physica
    finally:
        sys.setdlopenflags(dlflags)
    return physica

physica = importPhysica()
CXXObj = physica.CXXObj

del sys
del os

from . import Core
from .Core import *
