#!/usr/bin/env python3
"""Derive the ACTUAL third-party imports of the benchmark's Python tools and
try importing each one.

Iteration 004 lost a second array because the install list was written from
memory - numpy/scipy/pandas/torch/absl-py - while BEATRICE also imports
imageio, seaborn and tqdm, and FINEMAP-inf imports bgzip. Every fit failed with
ModuleNotFoundError after 33 minutes of queueing and a full re-run.

The source is on disk, so there is no reason to guess. This walks it with ast,
drops the standard library via sys.stdlib_module_names, drops modules that are
local to the tool tree itself, and reports what is genuinely missing BY NAME.

  py_deps_check.py <dir> [<dir> ...]

Exit 0 if every discovered import resolves, 1 otherwise.
"""
import ast
import importlib.util
import pathlib
import sys
import warnings

# Third-party sources routinely contain invalid escape sequences that emit
# SyntaxWarning when parsed. They are irrelevant here and only obscure the
# MISSING list, which is the one thing a reader of this output needs to see.
warnings.filterwarnings("ignore", category=SyntaxWarning)

# import-name -> pip-name, only where they differ.
PIP_NAME = {
    "absl": "absl-py", "sklearn": "scikit-learn", "cv2": "opencv-python",
    "PIL": "pillow", "yaml": "pyyaml", "skimage": "scikit-image",
    "Bio": "biopython", "dateutil": "python-dateutil",
}

roots = [pathlib.Path(a) for a in sys.argv[1:]]
roots = [r for r in roots if r.exists()]
if not roots:
    print("no readable tool directories given", file=sys.stderr)
    sys.exit(1)

found, local = set(), set()
for root in roots:
    for p in root.rglob("*.py"):
        local.add(p.stem)
        local.add(p.parent.name)
        try:
            tree = ast.parse(p.read_text(errors="ignore"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    found.add(a.name.split(".")[0])
            elif isinstance(node, ast.ImportFrom):
                # level > 0 is a relative import: local by definition.
                if node.level == 0 and node.module:
                    found.add(node.module.split(".")[0])

third_party = sorted(found - set(sys.stdlib_module_names) - local)

missing = []
for m in third_party:
    if importlib.util.find_spec(m) is None:
        missing.append(m)

print(f"third-party imports discovered: {len(third_party)}")
print("  " + " ".join(third_party))
if missing:
    print(f"\nMISSING ({len(missing)}):")
    for m in missing:
        print(f"  {m:<20} pip install {PIP_NAME.get(m, m)}")
    print("\npip install " + " ".join(PIP_NAME.get(m, m) for m in missing))
    sys.exit(1)
print("\nall discovered imports resolve")
