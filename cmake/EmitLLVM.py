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
import shutil
import subprocess
import sys
from pathlib import Path

def compile(device_only: bool):
    cuda_flag = "--offload-device-only" if device_only else "--offload-host-only"

    subprocess.run(["cmake", f"-DCMAKE_CUDA_FLAGS={cuda_flag}", ".."], stdout=subprocess.DEVNULL, check=True)
    subprocess.run(["cmake", "--build", ".", "--target=Benchmark"], check=True)

def collect(llvm_dir: Path, device_only: bool, arch: str):
    benchmark_dir = Path(".") / "benchmark"
    for file in benchmark_dir.rglob("*.o"):
        if "Dispatch" in str(file) or "MKL" in str(file):
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
    def get_arch(arch: str) -> str:
        arch = arch.strip()
        if arch == "native":
            def native() -> str:
                try:
                    result = subprocess.run(["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader"],
                                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
                except (OSError, subprocess.CalledProcessError):
                    return "native"

                result = result.stdout.split()
                return result[0].replace(".", "") if result else "native"
            arch = native()

        import re
        if arch == "" or not re.fullmatch(r"[0-9]+", arch):
            raise ValueError("Bad arch: {}".format(arch))
        return arch

    arch = get_arch(sys.argv[1]) if len(sys.argv) > 1 else get_arch("native")

    llvm_dir = Path(".") / "llvm"
    llvm_dir.mkdir(parents=True, exist_ok=True)

    def clean(llvm_dir: Path):
        for file in llvm_dir.glob("*.ll"):
            file.unlink()

    clean(llvm_dir)

    print("Collecting LLVM IR... ", end=" ")
    if arch == "":
        compile(False)
        collect(False, arch)
    else:
        for device_only in (False, True):
            compile(device_only)
            collect(llvm_dir, device_only, arch)
    print("Complete")
