"""Plan objectives, feasibility, Pareto front and hypervolume (stdlib only)."""
import random

OBJECTIVE_LABELS = {
    "coverageLoss": "Target coverage loss, nominal (% of Rx below D95)",
    "letBurden": "OAR high-LET burden",
    "robustLoss": "Target coverage loss, worst case (% of Rx below D95)",
}

LET_METRICS = {
    "dirtyDoseD2": ("dirtyDoseD2", "Dirty dose D2% summed over LET OARs (Gy)"),
    "VDLpct": ("VDLpct", "Volume with high dose and high LETd, summed over LET OARs (%)"),
    "LETdP98": ("LETdP98", "LETd near-max (98th percentile), max over LET OARs (keV/um)"),
}


def _as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def objectives(run, case, settings):
    """Objective vector (all minimised) and feasibility notes for one run."""
    m = run["metrics"]
    rx = case["prescription"]["totalDoseGy"]
    t = _as_list(m["targets"])[0]
    obj = {"coverageLoss": 100.0 * (1 - t["D95"]["nom"] / t["doseGy"])}
    key, _ = LET_METRICS[settings["letMetric"]]
    lets = [o for o in _as_list(m["oars"]) if o.get("let")]
    if key == "LETdP98":
        obj["letBurden"] = max([o[key]["nom"] for o in lets] or [0.0])
    else:
        obj["letBurden"] = sum(o[key]["nom"] for o in lets)
    if m.get("numOfScenarios", 1) > 1:
        obj["robustLoss"] = 100.0 * (1 - t["D95"]["wc"] / t["doseGy"])
    issues = []
    hard = {o["name"]: o for o in _as_list(case["structures"].get("oars")) if o.get("hardLimit")}
    for o in _as_list(m["oars"]):
        if o["name"] in hard and o.get("maxDoseGy") is not None:
            lim = o["maxDoseGy"] * (1 + settings["oarLimitTolerance"])
            if o["D2"]["nom"] > lim:
                issues.append("{} D2% {:.1f} Gy > limit {:.1f} Gy".format(o["name"], o["D2"]["nom"], o["maxDoseGy"]))
    hs = settings.get("maxHotspotFraction")
    if hs and t["D2"]["nom"] > hs * t["doseGy"]:
        issues.append("target D2% {:.1f} Gy > {:.0f}% of Rx".format(t["D2"]["nom"], 100 * hs))
    return obj, issues


def dominates(a, b, keys):
    return all(a[k] <= b[k] for k in keys) and any(a[k] < b[k] for k in keys)


def pareto_front(points, keys):
    """points: list of (id, objdict). Returns ids on the front."""
    front = []
    for i, (pid, a) in enumerate(points):
        if not any(dominates(b, a, keys) for j, (_, b) in enumerate(points) if j != i):
            front.append(pid)
    return front


def hypervolume(points, keys, bounds, samples=20000, seed=7):
    """Monte-Carlo hypervolume in normalised objective space.
    bounds: {key: (lo, hi)} fixed per search so values are comparable."""
    if not points:
        return 0.0
    norm = []
    for _, a in points:
        norm.append([(a[k] - bounds[k][0]) / max(1e-9, bounds[k][1] - bounds[k][0]) for k in keys])
    rng = random.Random(seed)
    ref = 1.1
    hit = 0
    for _ in range(samples):
        x = [rng.uniform(-0.1, ref) for _ in keys]
        if any(all(p[d] <= x[d] for d in range(len(keys))) for p in norm):
            hit += 1
    return hit / samples * (ref + 0.1) ** len(keys)
