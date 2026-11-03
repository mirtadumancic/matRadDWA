"""Diagnostic rules: the agent's judgement between rounds.

Each rule reads the history (recipes + objectives only), may change the
search space, and explains what it saw and did in plain words. These rules
encode what planners otherwise find by hand, e.g. LET objectives eroding
target coverage, or a LET objective that cannot act because no beam ends
near the OAR.
"""
import math
from .search import combo, penalty_range, LET_MODES


def _ok(h):
    return [r for r in h if r["status"] == "ok"]


def _anchor(history, bs):
    for r in _ok(history):
        p = r["params"]
        if (p["beamSet"] == bs and p["mode"] == "none" and p["robust"] == "none"
                and p["coverageMode"] == "penalty"):
            return r
    return None


def apply_rules(state, settings, front_ids):
    space, history = state["space"], state["history"]
    notes = []
    tol = settings["coverageErosionPct"]

    for bs in list(space["beamSets"]):
        anc = _anchor(history, bs)
        if anc is None:
            continue
        base_cov = anc["obj"]["coverageLoss"]
        base_let = anc["obj"]["letBurden"]
        for mode in LET_MODES:
            c = combo(bs, mode)
            if c in space["disabled"] or mode not in space["modes"]:
                continue
            runs = [r for r in _ok(history) if r["params"]["beamSet"] == bs and r["params"]["mode"] == mode
                    and r["params"]["robust"] == "none" and r["params"]["coverageMode"] == "penalty"]
            if not runs:
                continue

            # R1 coverage erosion: cap the penalty and allow a coverage constraint
            bad = [r for r in runs if r["obj"]["coverageLoss"] > base_cov + tol]
            good = [r for r in runs if r["obj"]["coverageLoss"] <= base_cov + tol]
            if bad and "constraint" not in space["coverageModesBy"].get(c, ["penalty"]):
                lo, hi = penalty_range(space, bs, mode)
                first_bad = min(r["params"]["penalty"] for r in bad)
                cap = 3 * max(r["params"]["penalty"] for r in good) if good else lo * 3
                cap = max(lo * 1.5, min(hi, first_bad, cap))
                space["penaltyRangeBy"][c] = [lo, cap]
                space["coverageModesBy"][c] = ["penalty", "constraint"]
                worst = max(bad, key=lambda r: r["obj"]["coverageLoss"])
                notes.append(dict(rule="coverage-erosion", text=(
                    "Beam set {bs}, {m}: at penalty {p:.3g} the target lost {d:.1f} more points of coverage than the "
                    "dose-only plan. Capped the LET penalty at {cap:.3g} and added plans that keep coverage as a hard "
                    "DVH constraint (V{f:.0f}%Rx >= 95%) instead of a weighted objective.").format(
                        bs=bs, m=mode, p=worst["params"]["penalty"], d=worst["obj"]["coverageLoss"] - base_cov,
                        cap=cap, f=95)))

            # R2 LET objective inactive
            pens = [r["params"]["penalty"] for r in runs]
            if len(runs) >= 3 and max(pens) / max(min(pens), 1e-9) >= 10:
                lets = [r["obj"]["letBurden"] for r in runs]
                change = (max(lets) - min(lets)) / max(abs(base_let), 1e-6)
                if base_let < settings["letNegligible"]:
                    space["disabled"].append(c)
                    notes.append(dict(rule="let-unneeded", text=(
                        "Beam set {bs}: the dose-only plan already has a negligible OAR LET burden ({v:.2f}); "
                        "stopped exploring {m} for this beam set.").format(bs=bs, v=base_let, m=mode)))
                elif change < 0.10:
                    if c not in space["expanded"]:
                        lo, hi = penalty_range(space, bs, mode)
                        space["penaltyRangeBy"][c] = [lo, hi * 10]
                        space["expanded"].append(c)
                        notes.append(dict(rule="let-inactive", text=(
                            "Beam set {bs}, {m}: a {r:.0f}x range of penalties changed the OAR LET burden by only "
                            "{ch:.0%}. Either the beams' distal edges do not end in the OAR or the penalty is too "
                            "weak; widened the penalty range 10x once.").format(
                                bs=bs, m=mode, r=max(pens) / min(pens), ch=change)))
                    else:
                        space["disabled"].append(c)
                        notes.append(dict(rule="let-inactive", text=(
                            "Beam set {bs}, {m}: still no LET response after widening the penalty range; dropped "
                            "this combination. A different beam arrangement is the likelier fix.").format(bs=bs, m=mode)))

    # R3 compare LET formulations on the front
    by_id = {r["runId"]: r for r in history}
    if all(m in space["modes"] for m in LET_MODES):
        counts = {m: sum(1 for i in front_ids if by_id[i]["params"]["mode"] == m) for m in LET_MODES}
        tried = {m: sum(1 for r in _ok(history) if r["params"]["mode"] == m) for m in LET_MODES}
        for m, other in (("LETxDose", "dirtyDose"), ("dirtyDose", "LETxDose")):
            if tried[m] >= settings["minRunsToJudge"] and tried[other] >= settings["minRunsToJudge"] \
                    and counts[m] == 0 and counts[other] >= 2:
                space["modes"].remove(m)
                notes.append(dict(rule="formulation", text=(
                    "No {m} plan is on the Pareto front after {n} tries, while {o} has {k}. Dropped {m} and "
                    "focused on {o}.{extra}").format(m=m, n=tried[m], o=other, k=counts[other],
                    extra=(" This matches Seckler et al. 2026, who found dirty dose gives more explicit LET control "
                           "than LET x dose." if other == "dirtyDose" else ""))))
                break

    # R4 robustness
    if space["robustAvailable"] and "COWC" not in space["robust"]:
        gaps = [by_id[i]["obj"].get("robustLoss", 0) - by_id[i]["obj"]["coverageLoss"] for i in front_ids]
        if gaps and sorted(gaps)[len(gaps) // 2] > settings["robustGapPct"]:
            space["robust"].append("COWC")
            notes.append(dict(rule="robustness", text=(
                "Front plans lose a median {g:.1f} more points of coverage in the worst-case scenario than "
                "nominally. Added robust (composite worst-case) plans to the search.").format(
                    g=sorted(gaps)[len(gaps) // 2])))

    # R5 beam sets that never reach the front
    if len(space["beamSets"]) > 1:
        for bs in list(space["beamSets"]):
            n_bs = sum(1 for r in _ok(history) if r["params"]["beamSet"] == bs)
            on_front = sum(1 for i in front_ids if by_id[i]["params"]["beamSet"] == bs)
            if n_bs >= settings["minRunsToJudge"] + 2 and on_front == 0 and len(space["beamSets"]) > 1:
                space["beamSets"].remove(bs)
                notes.append(dict(rule="beam-set", text=(
                    "Beam set {bs}: none of its {n} plans is on the Pareto front; stopped spending runs on it.")
                    .format(bs=bs, n=n_bs)))
    return notes
