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
from enum import Enum
from phypy import CXXObj

class ScalarOption(Enum):
    Float32 = 0
    Float64 = 1
    FloatMP = 2

class Real(CXXObj):
    __pDecl = None

    def __init__(self, Option: ScalarOption = ScalarOption.Float64) -> None:
        super().__init__(Real.__pDecl, Option.value)
        super().construct()
        self.__Option = Option

    def __float__(self):
        rtnTyName = str(self.Option).lower()
        return float(self.call(rtnTyName, 'toMachine'))

    def __repr__(self):
        return repr(float(self))

    @property
    def Option(self):
        return self.__Option

    @staticmethod
    def include() -> None:
        if (Real.__pDecl is None):
            from phypy import physica
            Real.__pDecl = physica.include('Core/Scalar/Real.h')

Real.include()

__all__ = [ScalarOption.__name__, Real.__name__]
