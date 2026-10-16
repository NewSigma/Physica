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
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

_SOURCE_DIR = Path(__file__).resolve().parent.parent
_BUILD_DIR = _SOURCE_DIR / "build"
_IR_DIR = _BUILD_DIR / "emit-llvm"

_SEED_ALLOW = frozenset({
    "CMAKE_BUILD_TYPE",
    "CMAKE_INSTALL_PREFIX",
    "CMAKE_PREFIX_PATH",
    "CMAKE_CUDA_ARCHITECTURES",
    "CMAKE_LINKER_TYPE",
})
_COMPILERS = ("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "CMAKE_CUDA_COMPILER")
_CACHE_ENTRY = re.compile(r"^(?P<name>[^/#:\n][^:\n]*):(?P<type>[^=\n]*)=(?P<value>.*)$")


def read_cache(directory: Path) -> tuple[dict[str, str], dict[str, str]]:
    assert directory in (_BUILD_DIR, _IR_DIR), directory
    cache = directory / "CMakeCache.txt"
    if not cache.is_file():
        raise FileNotFoundError(
            "No CMake cache at {}; configure that build tree first.".format(cache))

    options = {}
    compilers = {}
    for line in cache.read_text().splitlines():
        match = _CACHE_ENTRY.match(line)
        if match is None:
            continue
        name, type, value = match.group("name"), match.group("type"), match.group("value")
        if type in ("INTERNAL", "STATIC"):
            continue
        if name in _COMPILERS:
            compilers[name] = value
        elif name in _SEED_ALLOW or (not name.startswith("CMAKE_") and type != "UNINITIALIZED"):
            options[name] = value
    return options, compilers


def init_ir_dir(compilers: dict):
    if (_IR_DIR / "CMakeCache.txt").is_file():
        _, cached = read_cache(_IR_DIR)

        def matches(name: str, path: Path) -> bool:
            return name in cached and os.path.realpath(cached[name]) == os.path.realpath(path)

        if not all(matches(name, path) for name, path in compilers.items()):
            shutil.rmtree(_IR_DIR)

    _IR_DIR.mkdir(parents=True, exist_ok=True)


def cuda_root(options: dict) -> Path | None:
    for name in ("CUDA_CUDART", "CUDA_cudart_LIBRARY", "CUDA_cudart_static_LIBRARY"):
        library = options.get(name)
        if library is None or library.endswith("-NOTFOUND"):
            continue
        parent = Path(library).parent
        if parent.name in ("lib64", "lib"):
            return parent.parent
    return None


def compile(options: dict, compilers: dict, device_only: bool):
    cuda_flag = "--offload-device-only" if device_only else "--offload-host-only"
    configure = ["cmake", "-G", "Ninja"]
    if root := cuda_root(options):
        configure.append(f"-DCUDAToolkit_ROOT={root}")
    configure += [f"-D{name}={value}" for name, value in options.items()]
    configure += [f"-D{name}={value}" for name, value in compilers.items()]
    configure += [
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_BENCHMARKS=ON",
        "-DPHYSICA_EmitLLVM=ON",
        "-DPHYSICA_TRANSFORMS=OFF",
        f"-DCMAKE_CUDA_FLAGS={cuda_flag}",
    ]
    configure.append(str(_SOURCE_DIR.resolve()))
    subprocess.run(configure, cwd=_IR_DIR, stdout=subprocess.DEVNULL, check=True)

    subprocess.run(["cmake", "--build", ".", "--target=Benchmark"], cwd=_IR_DIR, check=True)


def collect(llvm_dir: Path, device_only: bool, arch: str):
    benchmark_dir = _IR_DIR / "benchmark"
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

    options, compilers = read_cache(_BUILD_DIR)  # Fail before creating anything if not configured.
    if args.clang is not None:
        clang = args.clang.resolve()
        clangxx = clang.parent / "clang++"
        if not clang.is_file():
            raise FileNotFoundError("clang not found: {}".format(clang))
        if not clangxx.is_file():
            raise FileNotFoundError("clang++ not found: {}".format(clangxx))
        compilers = {
            "CMAKE_C_COMPILER": clang,
            "CMAKE_CXX_COMPILER": clangxx,
            "CMAKE_CUDA_COMPILER": clangxx,
        }
    else:
        compilers = {name: Path(path) for name, path in compilers.items()}

    init_ir_dir(compilers)
    arch = get_arch("native")
    llvm_dir = _IR_DIR / "llvm"
    llvm_dir.mkdir(parents=True, exist_ok=True)

    def clean(llvm_dir: Path):
        for file in llvm_dir.glob("*.ll"):
            file.unlink()

    clean(llvm_dir)

    print("Collecting LLVM IR... ", end=" ")
    if arch == "":
        compile(options, compilers, False)
        collect(llvm_dir, False, "")
    else:
        for device_only in (False, True):
            compile(options, compilers, device_only)
            collect(llvm_dir, device_only, arch)
    print("Complete")


if __name__ == "__main__":
    sys.exit(main())
