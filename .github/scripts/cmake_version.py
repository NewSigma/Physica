# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
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
