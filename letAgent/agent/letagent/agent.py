"""The search agent: propose plans, run them, read metrics, apply rules, repeat."""
import datetime
import json
import os
import random

from . import pareto, privacy, report, rules, search

DEFAULT_SETTINGS = {
    "budget": 40,              # total plans
    "batch": 8,                # plans per round after the anchors
    "letMetric": "dirtyDoseD2",
    "letModes": ["dirtyDose", "LETxDose"],
    "oarLimitTolerance": 0.02,
    "maxHotspotFraction": None,
    "coverageErosionPct": 2.0,  # points of D95 coverage a LET plan may lose vs dose-only
    "letNegligible": 0.5,       # LET burden below which LET optimisation is pointless
    "robustGapPct": 3.0,
    "minRunsToJudge": 4,
    "stallRounds": 2,
    "stallGain": 0.01,
    "seed": 1,
}


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def resolve_case(case_file, matrad_root):
    """Copy the case file into its work folder with an absolute workDir so
    MATLAB and Python agree on where everything lives."""
    with open(case_file) as f:
        case = json.load(f)
    privacy.check_alias(case["alias"])
    # Refuse patient cases that would read or write inside the repository or an
    # AI-connected folder, before anything is written anywhere.
    privacy.check_data_location(case, case_file, matrad_root)
    work = case.get("workDir") or os.path.join(matrad_root, "letAgent_runs")
    if not os.path.isabs(work):
        work = os.path.join(os.path.dirname(os.path.abspath(case_file)), work)
    case["workDir"] = os.path.abspath(work)
    case_dir = os.path.join(case["workDir"], case["alias"])
    os.makedirs(case_dir, exist_ok=True)
    resolved = os.path.join(case_dir, "case_resolved.json")
    with open(resolved, "w") as f:
        json.dump(case, f, indent=2)
    return case, case_dir, resolved


def allowed_strings(case):
    s = {"ok", "failed", "none", "dirtyDose", "LETxDose", "COWC", "STOCH", "VWWC", "penalty", "constraint", "1", "0"}
    for t in case["structures"]["targets"]:
        s.add(t["name"])
    for o in case["structures"].get("oars", []):
        s.add(o["name"])
    for b in case["beamSets"]:
        s.add(b["id"])
    return s


class Agent:
    def __init__(self, case_file, engine_factory, matrad_root, settings=None):
        self.settings = dict(DEFAULT_SETTINGS, **(settings or {}))
        self.case, self.case_dir, self.case_file = resolve_case(case_file, matrad_root)
        self.engine = engine_factory(self.case_dir)
        self.state_file = os.path.join(self.case_dir, "agent_state.json")
        self.rng = random.Random(self.settings["seed"])
        self.allowed = allowed_strings(self.case)

    # ---------------------------------------------------------------- state
    def _load_state(self, summary):
        if os.path.exists(self.state_file):
            with open(self.state_file) as f:
                return json.load(f)
        return {"alias": self.case["alias"], "space": search.init_space(summary, self.settings),
                "history": [], "log": [], "hv": [], "round": 0, "summary": summary}

    def _save(self, state):
        with open(self.state_file, "w") as f:
            json.dump(state, f, indent=1)

    def _log(self, state, rule, text):
        entry = {"round": state["round"], "time": _now(), "rule": rule, "text": privacy.scrub_text(text, 600)}
        state["log"].append(entry)
        print("[agent r{}] {}: {}".format(state["round"], rule, entry["text"]))

    # ----------------------------------------------------------------- loop
    def run(self):
        summary = self.engine.prepare_case(self.case_file)
        state = self._load_state(summary)
        if not state["history"]:
            self._log(state, "start", "Case {} with beam sets {}, {} robustness scenario(s). Budget {} plans."
                      .format(self.case["alias"], ", ".join(state["space"]["beamSets"]),
                              summary.get("numOfScenarios", 1), self.settings["budget"]))
        seen = {search.key(r["params"]) for r in state["history"]}
        stall = 0
        while len(state["history"]) < self.settings["budget"]:
            state["round"] += 1
            left = self.settings["budget"] - len(state["history"])
            if not state["history"]:
                props = search.dedupe(search.anchors(state["space"]), seen)
                why = "anchors: a dose-only plan and one plan per LET formulation for each beam set"
            else:
                fp, fo = self._front_params(state)
                n = min(self.settings["batch"], left)
                n_ex = max(1, n // 2)
                props = search.dedupe(search.exploit(state["space"], fp, fo, n - n_ex, self.rng), seen)
                props += search.dedupe(search.explore(state["space"], n - len(props), self.rng), seen)
                why = "{} refinements near the front, {} exploratory".format(
                    sum(p["origin"] != "explore" for p in props), sum(p["origin"] == "explore" for p in props))
            props = props[:left]
            if not props:
                self._log(state, "stop", "No new plans to propose; the search space is exhausted.")
                break
            self._log(state, "propose", "Round {}: {} plans ({}).".format(state["round"], len(props), why))
            runs = []
            for p in props:
                rid = "r{:03d}".format(len(state["history"]) + len(runs) + 1)
                runs.append({"runId": rid, "beamSet": p["beamSet"], "recipe": search.recipe(p), "_p": p})
            params = {r["runId"]: r.pop("_p") for r in runs}
            results, secs = self.engine.run_batch(self.case_file, runs, "round{:02d}".format(state["round"]))
            for res in results:
                res = privacy.scrub_run(res, self.allowed)
                entry = {"runId": res["runId"], "round": state["round"], "params": params[res["runId"]],
                         "status": res.get("status", "failed"), "seconds": res.get("seconds")}
                if entry["status"] == "ok":
                    entry["metrics"] = res["metrics"]
                    entry["obj"], entry["issues"] = pareto.objectives(res, self.case, self.settings)
                    entry["notes"] = res.get("notes") or []
                else:
                    entry["error"] = privacy.scrub_text(res.get("error"))
                    self._log(state, "run-failed", "{}: {}".format(res["runId"], entry["error"]))
                state["history"].append(entry)
            front = self._front_ids(state)
            hv = self._hv(state, front)
            state["hv"].append(hv)
            self._log(state, "result", "Round {} took {:.0f} min; {} plans on the front; hypervolume {:.3f}."
                      .format(state["round"], secs / 60, len(front), hv))
            for note in rules.apply_rules(state, self.settings, front):
                self._log(state, note["rule"], note["text"])
            self._save(state)
            report.write(self.case_dir, self.case, state, front, self.settings)
            if len(state["hv"]) >= 3:
                prev = state["hv"][-2]
                gain = (hv - prev) / max(prev, 1e-9)
                stall = stall + 1 if gain < self.settings["stallGain"] else 0
                if stall >= self.settings["stallRounds"]:
                    self._log(state, "stop", "Front stopped improving (<{:.0%} hypervolume gain for {} rounds)."
                              .format(self.settings["stallGain"], stall))
                    break
        front = self._front_ids(state)
        self._log(state, "done", "{} plans evaluated; {} on the Pareto front. Report: report.html"
                  .format(len(state["history"]), len(front)))
        self._save(state)
        report.write(self.case_dir, self.case, state, front, self.settings)
        return state

    # --------------------------------------------------------------- helpers
    def _feasible(self, state):
        ok = [r for r in state["history"] if r["status"] == "ok"]
        feas = [r for r in ok if not r["issues"]]
        return feas if feas else ok

    def _keys(self, state):
        ok = [r for r in state["history"] if r["status"] == "ok"]
        return [k for k in ("coverageLoss", "letBurden", "robustLoss") if ok and k in ok[0]["obj"]]

    def _front_ids(self, state):
        pts = [(r["runId"], r["obj"]) for r in self._feasible(state)]
        return pareto.pareto_front(pts, self._keys(state))

    def _front_params(self, state):
        by = {r["runId"]: r for r in state["history"]}
        ids = self._front_ids(state)
        return [by[i]["params"] for i in ids], [by[i]["obj"] for i in ids]

    def _hv(self, state, front):
        keys = self._keys(state)
        ok = [r for r in state["history"] if r["status"] == "ok"]
        if "bounds" not in state:
            # fixed after the anchor round so later rounds are comparable
            state["bounds"] = {k: (min(r["obj"][k] for r in ok), max(r["obj"][k] for r in ok) or 1.0) for k in keys}
            for k in keys:
                lo, hi = state["bounds"][k]
                if hi - lo < 1e-6:
                    state["bounds"][k] = (lo, lo + 1.0)
        by = {r["runId"]: r for r in state["history"]}
        return pareto.hypervolume([(i, by[i]["obj"]) for i in front], keys, state["bounds"])
