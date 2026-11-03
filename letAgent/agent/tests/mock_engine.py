#!/usr/bin/env python3
"""Stand-in for MATLAB/Octave: mimics letAgent_prepareCase / letAgent_runBatch with a toy dose/LET trade-off model."""
import sys, re, json, os, math, random
code = sys.argv[-1]
m = re.search(r"letAgent_prepareCase\('([^']+)'\)", code)
if m:
    case = json.load(open(m.group(1)))
    d = os.path.join(case["workDir"], case["alias"])
    nS = 9 if case["robustness"]["scenarios"] == "wcScen" else 1
    json.dump({"alias": case["alias"], "numOfScenarios": nS,
               "beamSets": [{"id": b["id"], "gantryAngles": b["gantryAngles"], "numOfBixels": 1000, "numOfScenarios": nS} for b in case["beamSets"]]},
              open(os.path.join(d, "case_summary.json"), "w"))
    sys.exit(0)
m = re.search(r"letAgent_runBatch\('([^']+)','([^']+)'\)", code)
batch = json.load(open(m.group(1))); case = json.load(open(batch["caseFile"]))
nS = 9 if case["robustness"]["scenarios"] == "wcScen" else 1
rx = case["prescription"]["totalDoseGy"]
# toy model: per beam set baseline dirty dose in brainstem, LET responsiveness, coverage cost
BS = {"lateralPair": (6.0, 1.0, 1.0), "awayFromBS": (0.3, 0.2, 0.5), "threeField": (4.0, 0.9, 0.8)}
out = []
for r in batch["runs"]:
    base, resp, cost = BS.get(r["beamSet"], (3, 0.8, 1))
    rec = r["recipe"]; mode = rec["let"]["mode"]; p = rec["let"].get("penalty", 0)
    oar = rec.get("oarDosePenaltyScale", 1)
    covloss, dd = 1.0 + 0.3 * math.log2(max(oar, 1e-3)) ** 2 * 0.3, base / (oar ** 0.15)
    if mode != "none":
        lo = {"dirtyDose": 10, "LETxDose": 1}[mode]
        x = max(0, math.log10(p / lo))
        eff = resp * (0.85 if mode == "dirtyDose" else 0.6)
        dd *= math.exp(-eff * x)
        c = cost * 0.25 * x ** 2.2 * (1.0 if mode == "dirtyDose" else 1.6)
        if rec["coverage"]["mode"] == "constraint":
            c = min(c, 0.8); dd *= 1.15
        covloss += c
    rob = rec.get("robust", "none")
    wcextra = 2.0 if rob == "COWC" else 5.5
    if rob == "COWC": covloss += 0.7; dd *= 0.9
    D95 = rx * (1 - covloss / 100); D95wc = rx * (1 - (covloss + wcextra) / 100)
    rnd = random.Random(r["runId"])
    dd *= 1 + rnd.uniform(-0.03, 0.03)
    out.append({"runId": r["runId"], "beamSet": r["beamSet"], "recipe": rec, "status": "ok", "notes": [], "seconds": 1.0,
                "metrics": {"numOfScenarios": nS, "prescriptionGy": rx,
                            "targets": [{"name": "PTV", "doseGy": rx, "D95": {"nom": D95, "wc": D95wc}, "D2": {"nom": rx * 1.05, "wc": rx * 1.08}, "V95pct": {"nom": 100 - covloss * 2, "wc": 90.0}}],
                            "oars": [{"name": "Brainstem", "maxDoseGy": 54, "let": True, "D2": {"nom": 50.0, "wc": 53.0}, "Dmean": {"nom": 20.0, "wc": 22.0},
                                      "dirtyDoseD2": {"nom": dd, "wc": dd * 1.3}, "VDLpct": {"nom": dd * 2, "wc": dd * 2.5}, "LETdP98": {"nom": 2 + dd / 3, "wc": 2.2 + dd / 3}},
                                     {"name": "BODY", "maxDoseGy": 30, "let": False, "D2": {"nom": 40.0, "wc": 42.0}, "Dmean": {"nom": 5.0, "wc": 5.5},
                                      "dirtyDoseD2": {"nom": 1, "wc": 1}, "VDLpct": {"nom": 1, "wc": 1}, "LETdP98": {"nom": 3, "wc": 3}}]}})
json.dump({"runs": out}, open(m.group(2), "w"))
