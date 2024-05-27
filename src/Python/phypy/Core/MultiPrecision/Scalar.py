"""
Copyright 2024 WeiBo He.

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
from enum import Enum
from phypy import CXXObj

class ScalarOption(Enum):
    Float = 0
    Double = 1
    MultiPrecision = 2

class Scalar(CXXObj):
    __pDecl = None

    def __init__(self, Option: ScalarOption) -> None:
        self.__Option = Option

    def __float__(self):
        return float(self.call(float, 'getTrivial'))

    def __repr__(self):
        return repr(float(self))

    @property
    def Option(self):
        return self.__Option

    @staticmethod
    def getPTU():
        pDecl = Scalar.__pDecl
        if (pDecl is None):
            from phypy import physica
            pDecl = physica.include('Core/MultiPrecision/Scalar.h')
        return pDecl
