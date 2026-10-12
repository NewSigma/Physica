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

    Keep it simple stupid:
    The build directory is `build/` at the repository root, the generator is Ninja, and the CUDA arch is native.
    The only knob is the clang path; `clang++` is taken from the same directory.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

_SOURCE_DIR = Path(__file__).resolve().parent.parent
_BUILD_DIR = _SOURCE_DIR / "build"


def compile(build_dir: Path, source_dir: Path, clang: Path | None, device_only: bool):
    cuda_flag = "--offload-device-only" if device_only else "--offload-host-only"

    configure = [
        "cmake",
        "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DPHYSICA_EmitLLVM=ON",
        "-DPHYSICA_TRANSFORMS=OFF",
        f"-DCMAKE_CUDA_FLAGS={cuda_flag}",
    ]
    if clang is not None:
        clangxx = clang.parent / "clang++"
        configure += [
            f"-DCMAKE_C_COMPILER={clang}",
            f"-DCMAKE_CXX_COMPILER={clangxx}",
            f"-DCMAKE_CUDA_COMPILER={clangxx}",
        ]
    configure.append(str(source_dir.resolve()))
    subprocess.run(configure, cwd=build_dir, stdout=subprocess.DEVNULL, check=True)
    subprocess.run(["cmake", "--build", ".", "--target=Benchmark"],
                   cwd=build_dir, check=True)


def collect(llvm_dir: Path, build_dir: Path, device_only: bool, arch: str):
    benchmark_dir = build_dir / "benchmark"
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

    if arch == "" or not re.fullmatch(r"[0-9]+", arch):
        raise ValueError("Bad arch: {}".format(arch))
    return arch


def main():
    parser = argparse.ArgumentParser(description="Generate LLVM IR for benchmarks")
    parser.add_argument("--clang", type=Path, default=None,
                        help="clang path; clang++ is taken from the same directory")
    args = parser.parse_args()

    if args.clang is not None:
        clang = args.clang.resolve()
        clangxx = clang.parent / "clang++"
        if not clang.is_file():
            raise FileNotFoundError("clang not found: {}".format(clang))
        if not clangxx.is_file():
            raise FileNotFoundError("clang++ not found: {}".format(clangxx))
    else:
        clang = None

    arch = get_arch("native")
    llvm_dir = _BUILD_DIR / "llvm"
    llvm_dir.mkdir(parents=True, exist_ok=True)

    def clean(llvm_dir: Path):
        for file in llvm_dir.glob("*.ll"):
            file.unlink()

    clean(llvm_dir)

    print("Collecting LLVM IR... ", end=" ")
    if arch == "":
        compile(_BUILD_DIR, _SOURCE_DIR, clang, False)
        collect(llvm_dir, _BUILD_DIR, False, "")
    else:
        for device_only in (False, True):
            compile(_BUILD_DIR, _SOURCE_DIR, clang, device_only)
            collect(llvm_dir, _BUILD_DIR, device_only, arch)
    print("Complete")


if __name__ == "__main__":
    sys.exit(main())
