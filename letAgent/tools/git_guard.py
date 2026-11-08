#!/usr/bin/env python3
"""Pre-commit guard: stop patient data and run outputs from being committed.

Blocks staged files that are
  * DICOM (.dcm / .ima, or files starting with the DICOM 'DICM' preamble),
  * letAgent run outputs (anything under a letAgent_runs* folder),
  * .mat / .raw / .mhd / .nii files outside matRad's machine and phantom folders,
  * letAgent case JSONs that are patient cases (dataClass 'patient', or a .mat source
    without dataClass).

Install once per clone:   python3 letAgent/tools/git_guard.py --install
(or, with the pre-commit framework, the 'letagent-guard' hook in .pre-commit-config.yaml)
"""
import json
import os
import re
import subprocess
import sys

ALLOWED_DATA_DIRS = ("matRad/basedata/", "matRad/phantoms/", "userdata/machines/", "test/",
                     "DWA_final_machine/", "IBA_data/")
DATA_EXT = (".mat", ".raw", ".mhd", ".nii", ".nii.gz", ".nrrd")


def staged():
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
                         capture_output=True, text=True, check=True).stdout
    return [l for l in out.splitlines() if l.strip()]


def is_dicom(path):
    if path.lower().endswith((".dcm", ".ima", ".dicom")):
        return True
    try:
        with open(path, "rb") as f:
            head = f.read(132)
        return len(head) == 132 and head[128:132] == b"DICM"
    except OSError:
        return False


def problem(path):
    p = path.replace("\\", "/")
    if re.search(r"(^|/)letAgent_runs[^/]*/", p):
        return "letAgent run output"
    if is_dicom(path):
        return "DICOM file"
    if p.lower().endswith(DATA_EXT) and not p.startswith(ALLOWED_DATA_DIRS):
        return "image/data file outside the machine and phantom folders"
    if p.startswith("letAgent/") and p.endswith(".json") and os.path.exists(path):
        try:
            with open(path) as f:
                c = json.load(f)
        except (ValueError, OSError):
            return None
        if isinstance(c, dict) and "alias" in c and "source" in c:
            dc = c.get("dataClass") or ("patient" if c.get("source", {}).get("type") == "matFile" else "synthetic")
            if dc == "patient" and not p.endswith("UHN_CNS_TEMPLATE.json"):
                return "patient case file"
    return None


def install():
    root = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                          check=True).stdout.strip()
    hook = os.path.join(root, ".git", "hooks", "pre-commit")
    if os.path.exists(hook):
        with open(hook) as f:
            if "git_guard.py" in f.read():
                print("letAgent guard already installed.")
                return 0
        print("A pre-commit hook already exists at {}; add this line to it:\n"
              "  python3 letAgent/tools/git_guard.py || exit 1".format(hook))
        return 1
    with open(hook, "w") as f:
        f.write("#!/bin/sh\npython3 letAgent/tools/git_guard.py || exit 1\n")
    os.chmod(hook, 0o755)
    print("Installed letAgent guard in", hook)
    return 0


def main():
    if "--install" in sys.argv:
        return install()
    bad = [(f, why) for f in staged() for why in [problem(f)] if why]
    if not bad:
        return 0
    print("letAgent guard: commit blocked. These staged files may hold patient data or run outputs:")
    for f, why in bad:
        print("  {}  ({})".format(f, why))
    print("Unstage them (git restore --staged <file>) and keep patient data on hospital storage.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
