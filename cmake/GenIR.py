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
    Generates LLVM IR for benchmarks. It is driven by CMake.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

def compile(device_only: bool):
    cuda_flag = ""
    if device_only:
        cuda_flag = "--offload-device-only"

    subprocess.run(["cmake", f"-DCMAKE_CUDA_FLAGS={cuda_flag}", ".."], stdout=subprocess.DEVNULL, check=True)
    subprocess.run(["cmake", "--build", ".", "--target=Benchmark"])

def collect(llvm_dir: Path, device_only: bool, arch: str):
    benchmark_dir = Path(".") / "benchmark"
    for file in benchmark_dir.rglob("*.o"):
        if "MKL" in str(file):
            continue

        if device_only:
            if (file.name.endswith(".cu.o")):
                shutil.move(file, llvm_dir / file.name.replace(".cu.o", "_cu_sm{}.ll".format(arch)))
            continue

        if file.name.endswith(".cpp.o"):
            shutil.move(file, llvm_dir / file.name.replace(".cpp.o", ".ll"))
        elif file.name.endswith(".cu.o"):
            shutil.move(file, llvm_dir / file.name.replace(".cu.o", "_cu.ll"))

if __name__ == "__main__":
    arch = ""
    if len(sys.argv) > 1:
        arch = sys.argv[1]
        if not arch.isnumeric():
            raise ValueError("Bad arch")

    llvm_dir = Path(".") / "llvm"
    llvm_dir.mkdir(parents=True, exist_ok=True)

    print("Collecting LLVM IR... ", end=" ")
    if arch == "":
        compile(False)
        collect(False, arch)
    else:
        for device_only in (False, True):
            compile(device_only)
            collect(llvm_dir, device_only, arch)
    print("Complete")
