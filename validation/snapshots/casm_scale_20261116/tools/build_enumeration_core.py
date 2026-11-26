"""Build the optional standard-library-only C++17 shared library in place."""
import argparse
from pathlib import Path
import platform
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="c++")
    args = parser.parse_args()
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    if sys.platform == "win32":
        parser.error("This build helper currently supports macOS/Linux")
    target = ROOT / "native" / ("libmof_enumeration" + suffix)
    command = [*shlex.split(args.compiler), "-std=c++17", "-O3", "-fPIC",
               "-dynamiclib" if sys.platform == "darwin" else "-shared",
               str(ROOT / "native" / "enumeration_core.cpp"), "-o", str(target)]
    subprocess.run(command, check=True)
    print(f"Built {target} ({platform.machine()})")
