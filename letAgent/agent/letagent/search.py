"""Search space and plan proposals (anchors, exploration, exploitation)."""
import math
import random

LET_MODES = ("dirtyDose", "LETxDose")
DEFAULT_PENALTY_RANGE = {"dirtyDose": [10.0, 3000.0], "LETxDose": [1.0, 300.0]}


def combo(bs, mode):
    return "{}|{}".format(bs, mode)


def init_space(summary, settings):
    beam_sets = [b["id"] for b in _as_list(summary["beamSets"])]
    return {
        "beamSets": beam_sets,
        "modes": ["none"] + list(settings.get("letModes", LET_MODES)),
        "penaltyRange": {m: list(DEFAULT_PENALTY_RANGE[m]) for m in LET_MODES},
        "penaltyRangeBy": {},
        "oarScaleRange": [0.25, 4.0],
        "robust": ["none"],
        "robustAvailable": summary.get("numOfScenarios", 1) > 1,
        "coverageModesBy": {},
        "disabled": [],
        "expanded": [],
    }


def _as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def penalty_range(space, bs, mode):
    return space["penaltyRangeBy"].get(combo(bs, mode), space["penaltyRange"][mode])


def coverage_modes(space, bs, mode):
    return space["coverageModesBy"].get(combo(bs, mode), ["penalty"])


def allowed_modes(space, bs):
    return [m for m in space["modes"] if combo(bs, m) not in space["disabled"]]


def recipe(p):
    r = {"let": {"mode": p["mode"]}, "oarDosePenaltyScale": round(p["oarScale"], 4),
         "robust": p["robust"], "coverage": {"mode": p["coverageMode"]}}
    if p["mode"] != "none":
        r["let"]["penalty"] = round(p["penalty"], 4)
    return r


def key(p):
    pen = round(math.log10(p["penalty"]), 1) if p["mode"] != "none" else None
    return (p["beamSet"], p["mode"], pen, round(math.log10(p["oarScale"]), 1), p["robust"], p["coverageMode"])


def anchors(space):
    out = []
    for bs in space["beamSets"]:
        out.append(dict(beamSet=bs, mode="none", penalty=0, oarScale=1.0, robust="none",
                        coverageMode="penalty", origin="anchor"))
        for m in allowed_modes(space, bs):
            if m == "none":
                continue
            lo, hi = penalty_range(space, bs, m)
            out.append(dict(beamSet=bs, mode=m, penalty=math.sqrt(lo * hi), oarScale=1.0,
                            robust="none", coverageMode="penalty", origin="anchor"))
    return out


def _logu(rng, lo, hi):
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def explore(space, n, rng):
    out = []
    for _ in range(n):
        bs = rng.choice(space["beamSets"])
        modes = allowed_modes(space, bs)
        let_modes = [m for m in modes if m != "none"]
        mode = rng.choice(let_modes) if let_modes and rng.random() > 0.15 else "none"
        lo, hi = penalty_range(space, bs, mode) if mode != "none" else (1, 1)
        out.append(dict(beamSet=bs, mode=mode, penalty=_logu(rng, lo, hi) if mode != "none" else 0,
                        oarScale=_logu(rng, *space["oarScaleRange"]), robust=rng.choice(space["robust"]),
                        coverageMode=rng.choice(coverage_modes(space, bs, mode)), origin="explore"))
    return out


def exploit(space, front_params, front_objs, n, rng):
    """Perturb front members and fill the largest gaps along the front."""
    out = []
    if not front_params:
        return out
    # gap filling: neighbours along the LET axis that share beam set and mode
    order = sorted(range(len(front_params)), key=lambda i: front_objs[i]["letBurden"])
    gaps = []
    for a, b in zip(order, order[1:]):
        pa, pb = front_params[a], front_params[b]
        if pa["beamSet"] == pb["beamSet"] and pa["mode"] == pb["mode"] and pa["mode"] != "none":
            d = sum(abs(front_objs[a][k] - front_objs[b][k]) for k in front_objs[a])
            gaps.append((d, pa, pb))
    for d, pa, pb in sorted(gaps, key=lambda g: -g[0])[: max(1, n // 3)]:
        q = dict(pa)
        q["penalty"] = math.sqrt(pa["penalty"] * pb["penalty"])
        q["oarScale"] = math.sqrt(pa["oarScale"] * pb["oarScale"])
        q["origin"] = "gap"
        out.append(q)
    i = 0
    while len(out) < n:
        p = dict(front_params[i % len(front_params)])
        i += 1
        if p["beamSet"] not in space["beamSets"] or combo(p["beamSet"], p["mode"]) in space["disabled"]:
            if i > 4 * len(front_params):
                break
            continue
        if p["mode"] != "none":
            lo, hi = penalty_range(space, p["beamSet"], p["mode"])
            p["penalty"] = min(hi, max(lo, p["penalty"] * math.exp(rng.gauss(0, 0.6))))
        lo, hi = space["oarScaleRange"]
        p["oarScale"] = min(hi, max(lo, p["oarScale"] * math.exp(rng.gauss(0, 0.3))))
        if len(space["robust"]) > 1 and rng.random() < 0.3:
            p["robust"] = rng.choice(space["robust"])
        cm = coverage_modes(space, p["beamSet"], p["mode"])
        if len(cm) > 1 and rng.random() < 0.3:
            p["coverageMode"] = rng.choice(cm)
        p["origin"] = "refine"
        out.append(p)
    return out


def dedupe(proposals, seen):
    out = []
    for p in proposals:
        k = key(p)
        if k not in seen:
            seen.add(k)
            out.append(p)
    return out
