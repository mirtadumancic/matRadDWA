"""Command line:

    python -m letagent search --case ../cases/TG119.json --engine matlab --budget 40 --workers 2
    python -m letagent prepare --case ../cases/TG119.json --engine matlab

Run from letAgent/agent, or put letAgent/agent on PYTHONPATH.
"""
import argparse
import os
import sys

from .agent import Agent, DEFAULT_SETTINGS, resolve_case
from .engine import Engine


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    default_root = os.path.abspath(os.path.join(here, "..", "..", ".."))
    ap = argparse.ArgumentParser(prog="letagent", description="LET-aware multi-criteria plan search on matRad")
    ap.add_argument("command", choices=["prepare", "search"])
    ap.add_argument("--case", required=True, help="case JSON file")
    ap.add_argument("--engine", choices=["matlab", "octave"], default="matlab")
    ap.add_argument("--engine-cmd", default=None, help="path to the matlab/octave executable if not on PATH")
    ap.add_argument("--matrad-root", default=default_root, help="matRad folder (default: this repository)")
    ap.add_argument("--workers", type=int, default=1, help="parallel MATLAB processes (each loads the dij into RAM)")
    ap.add_argument("--budget", type=int, default=DEFAULT_SETTINGS["budget"])
    ap.add_argument("--batch", type=int, default=DEFAULT_SETTINGS["batch"])
    ap.add_argument("--let-metric", choices=["dirtyDoseD2", "VDLpct", "LETdP98"], default=DEFAULT_SETTINGS["letMetric"])
    ap.add_argument("--seed", type=int, default=DEFAULT_SETTINGS["seed"])
    a = ap.parse_args(argv)

    def factory(case_dir):
        return Engine(a.engine, a.matrad_root, case_dir, matlab_cmd=a.engine_cmd, workers=a.workers)

    if a.command == "prepare":
        case, case_dir, resolved = resolve_case(a.case, a.matrad_root)
        s = factory(case_dir).prepare_case(resolved)
        print("Prepared {}: beam sets {}".format(case["alias"], [b["id"] for b in (s["beamSets"] if isinstance(s["beamSets"], list) else [s["beamSets"]])]))
        return 0
    settings = {"budget": a.budget, "batch": a.batch, "letMetric": a.let_metric, "seed": a.seed}
    Agent(a.case, factory, a.matrad_root, settings).run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
