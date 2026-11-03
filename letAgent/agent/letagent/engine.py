"""Runs the matRad side (MATLAB or Octave) as subprocesses.

Only JSON goes in and out: a batch of recipes in, metrics out. Large data
(CT, influence matrices, spot weights) stays in the case work directory and
is never read by the agent.
"""
import json
import os
import subprocess
import time


def _mpath(p):
    """Path as a MATLAB string literal (forward slashes work on Windows too)."""
    return os.path.abspath(p).replace("\\", "/").replace("'", "''")


class Engine:
    def __init__(self, kind, matrad_root, work_case_dir, matlab_cmd=None, workers=1):
        if kind not in ("matlab", "octave"):
            raise ValueError("engine must be 'matlab' or 'octave'")
        self.kind = kind
        self.matrad_root = os.path.abspath(matrad_root)
        self.letagent_m = os.path.join(self.matrad_root, "letAgent", "matlab")
        self.case_dir = work_case_dir
        self.log_dir = os.path.join(work_case_dir, "logs")
        os.makedirs(self.log_dir, exist_ok=True)
        self.cmd = matlab_cmd or ("matlab" if kind == "matlab" else "octave")
        self.workers = max(1, int(workers))

    def _code(self, call):
        return ("cd('{root}'); matRad_rc; addpath('{la}'); "
                "cfg = MatRad_Config.instance(); cfg.disableGUI = true; cfg.logLevel = 1; "
                "{call}").format(root=_mpath(self.matrad_root), la=_mpath(self.letagent_m), call=call)

    def _argv(self, code):
        if self.kind == "matlab":
            return [self.cmd, "-batch", code]
        return [self.cmd, "--no-gui", "--quiet", "--eval", code]

    def _start(self, call, log_name):
        log = open(os.path.join(self.log_dir, log_name), "w")
        p = subprocess.Popen(self._argv(self._code(call)), stdout=log, stderr=subprocess.STDOUT,
                             cwd=self.matrad_root)
        return p, log

    def prepare_case(self, case_file):
        p, log = self._start("letAgent_prepareCase('{}');".format(_mpath(case_file)), "prepare.log")
        rc = p.wait()
        log.close()
        if rc != 0:
            raise RuntimeError("Case preparation failed; see logs/prepare.log in the case folder.")
        with open(os.path.join(self.case_dir, "case_summary.json")) as f:
            return json.load(f)

    def run_batch(self, case_file, runs, tag):
        """Run recipes on up to `workers` parallel processes. Returns results
        in the order of `runs`. Runs of the same beam set go to the same
        worker so each worker loads each cache once."""
        runs = sorted(runs, key=lambda r: r["beamSet"])
        n = min(self.workers, len(runs))
        chunks = [runs[i * len(runs) // n:(i + 1) * len(runs) // n] for i in range(n)]
        procs = []
        for i, chunk in enumerate(chunks):
            bfile = os.path.join(self.case_dir, "batches", "{}_w{}.json".format(tag, i))
            ofile = bfile.replace(".json", "_results.json")
            os.makedirs(os.path.dirname(bfile), exist_ok=True)
            with open(bfile, "w") as f:
                json.dump({"caseFile": os.path.abspath(case_file), "runs": chunk}, f, indent=1)
            if os.path.exists(ofile):
                os.remove(ofile)
            p, log = self._start("letAgent_runBatch('{}','{}');".format(_mpath(bfile), _mpath(ofile)),
                                 "{}_w{}.log".format(tag, i))
            procs.append((p, log, ofile, chunk))
        t0 = time.time()
        results = {}
        for p, log, ofile, chunk in procs:
            p.wait()
            log.close()
            got = {}
            if os.path.exists(ofile):
                with open(ofile) as f:
                    data = json.load(f)
                rr = data.get("runs", [])
                if isinstance(rr, dict):
                    rr = [rr]
                got = {r["runId"]: r for r in rr}
            for r in chunk:
                results[r["runId"]] = got.get(r["runId"], dict(r, status="failed",
                                              error="worker stopped before this run finished"))
        return [results[r["runId"]] for r in runs], time.time() - t0
