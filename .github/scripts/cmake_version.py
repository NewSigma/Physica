"""
Copyright 2026 Weibo He.

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

"""
    Reports the minimum CMake version declared by `CMakeLists.txt`.

    CI installs exactly this version, so that every build exercises the lowest supported CMake.
    Keep it simple stupid: the version has a single source, the declaration in `CMakeLists.txt`.
"""
import re
from pathlib import Path

_SOURCE_DIR = Path(__file__).resolve().parents[2]
_DECLARATION = re.compile(r"^\s*cmake_minimum_required\s*\(\s*VERSION\s+(\d+\.\d+)", re.MULTILINE)


if __name__ == "__main__":
    cmake_lists = _SOURCE_DIR / "CMakeLists.txt"
    match = _DECLARATION.search(cmake_lists.read_text(encoding="utf-8"))
    if match is None:
        raise SystemExit("Cannot parse cmake_minimum_required from {}".format(cmake_lists))
    print(match.group(1))
