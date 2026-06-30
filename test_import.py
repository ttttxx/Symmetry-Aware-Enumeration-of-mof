"""
Quick import test for the Streamlit frontend dependencies.
Run with: .venv/bin/python test_import.py
"""
import sys
print(f"Python: {sys.executable}  ({sys.version.split()[0]})")

results = {}

def try_import(name, from_str=None):
    try:
        if from_str:
            mod = __import__(from_str, fromlist=[name])
            obj = getattr(mod, name)
        else:
            mod = __import__(name)
        results[name] = "OK"
    except Exception as e:
        results[name] = f"FAIL — {e}"

# Core scientific
try_import("numpy")
try_import("pymatgen")
try_import("Structure", from_str="pymatgen.core")

# 3-D visualisation
try_import("py3Dmol")
try_import("stmol")

# Streamlit
try_import("streamlit")

# Project modules (relative to this file's directory)
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try_import("symmetry")
try_import("enumerator")
try_import("verification")
try_import("assembler")

try:
    from main import run_enumeration
    results["main.run_enumeration"] = "OK"
except Exception as e:
    results["main.run_enumeration"] = f"FAIL — {e}"

# ── Report ────────────────────────────────────────────────────────────────────
print()
all_ok = True
for pkg, status in results.items():
    mark = "✓" if status == "OK" else "✗"
    print(f"  {mark}  {pkg:<30} {status}")
    if status != "OK":
        all_ok = False

print()
if all_ok:
    print("All imports OK — you can now run:  .venv/bin/streamlit run app.py")
else:
    print("Some imports FAILED — see above.")
    sys.exit(1)
